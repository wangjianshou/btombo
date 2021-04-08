import io
from tombo import tombo_helper as th
from tombo._preprocess import _prep_fast5_for_fastq # 获取single_fast5的read_id,创建注释slot 
from tombo._preprocess import _get_ann_queues
from concurrent.futures import ThreadPoolExecutor, as_completed

VERBOSE = False

_MAX_QUEUE_SIZE = 1000
_ITER_QUEUE_LIMIT = 1000
_PROC_UPDATE_INTERVAL = 100

_MAX_FASTQ_QUEUE_SIZE = 10000
_SEQ_SUMMARY_FN_FIELD = 'filename'
_SEQ_SUMMARY_ID_FIELD = 'read_id'

# warning messages for annotate with fastqs over multiple processes,
# requiring passing warning codes to only print warning once.
_WARN_ID_VAL = 'ids'
_WARN_IO_VAL = 'io'
_WARN_MISMATCH_VAL = 'mismatch'
_WARN_OVRWRT_VAL = 'overwrite'
_WARN_UNIQ_VAL = 'uniq'
_WARN_CODES = (_WARN_ID_VAL, _WARN_IO_VAL, _WARN_MISMATCH_VAL, _WARN_OVRWRT_VAL)
_WARN_CODES_PREP = (_WARN_OVRWRT_VAL, _WARN_UNIQ_VAL)
_WARN_PREFIX = '****** WARNING ****** '


#单线程版，弃用
def get_seq_records(fastq_fns):
    fastq_recs = {}
    for fastq_fn in fastq_fns:
        with io.open(fastq_fn) as fastq_fp:
            tmp = [i for i in fastq_fp]
        if len(tmp) % 4 != 0: break
        tmp = {tmp[4*i].split()[0].split('_')[0][1:]:tmp[4*i:4*i+4]
               for i in range(len(tmp)//4)}
        for k in tmp:
            if not (tmp[k][0].startswith('@') and tmp[k][2].startswith('+')):
                th.warning_message(
                        'Successfully parsed ' + unicode(n_recs) +
                        ' FASTQ records from ' + fastq_fn + ' before ' +
                        'encountering an invalid record. The rest of ' +
                        'this file will not be processed.')
        fastq_recs.update(tmp)
    return fastq_recs


# 多线程版
def get_seq_recs_concurrent(fastq_fns, n=6):
    #fastq_recs = {}
    exector = ThreadPoolExecutor(max_workers=n)
    all_tasks = [exector.submit(get_seq_worker, fastq_fn) for fastq_fn in fastq_fns]
    return [future.result() for future in as_completed(all_tasks)]
    #for future in as_completed(all_tasks):
    #    fastq_recs.update(future.result())
    #return fastq_recs


def get_seq_worker(fastq_fn):
    with io.open(fastq_fn) as fastq_fp:
        tmp = [i for i in fastq_fp]
    if len(tmp) % 4 != 0: 
        raise Exception('error in fastq file')
    tmp = {tmp[4*i].split()[0].split('_')[0][1:]:tmp[4*i:4*i+4]
                for i in range(len(tmp)//4)}
    for k in tmp:
        if not (tmp[k][0].startswith('@') and tmp[k][2].startswith('+')):
            th.warning_message(
                    'Successfully parsed ' + unicode(n_recs) +
                    ' FASTQ records from ' + fastq_fn + ' before ' +
                    'encountering an invalid record. The rest of ' +
                    'this file will not be processed.')
    return tmp



def _annotate_with_fastqs_worker(single_fast5_q, fastq_recs, fastq_slot,  
            prog_q, warn_q, bc_grp_name, bc_subgrp_name, overwrite): 
    been_warned = dict((warn_code, False) for warn_code in _WARN_CODES)
    num_recs_proc = 0
    while True:
        single_fast5 = single_fast5_q.get()
        if single_fast5 is None:
            break
        try:
            with h5py.File(single_fast5, 'r+') as fast5_data:
                #if not fq_slot_prepped:
                try:
                    file_parsed_id = _prep_fast5_for_fastq(
                        fast5_data, bc_grp_name, bc_subgrp_name, overwrite)
                except th.TomboError:
                    if not been_warned[_WARN_OVRWRT_VAL]:
                        been_warned[_WARN_OVRWRT_VAL] = True
                        warn_q.put(_WARN_OVRWRT_VAL)
                    continue
                if file_parsed_id not in fastq_recs:
                    if not been_warned[_WARN_MISMATCH_VAL]:
                        been_warned[_WARN_MISMATCH_VAL] = True
                        warn_q.put(_WARN_MISMATCH_VAL)
                    continue
                bc_slot = fast5_data[fastq_slot]
                # add sequence to fastq slot
                bc_slot.create_dataset(
                    'Fastq', data=''.join(fastq_rec),
                    dtype=h5py.special_dtype(vlen=unicode))
                num_recs_proc += 1
                if num_recs_proc % _PROC_UPDATE_INTERVAL == 0:
                    prog_q.put(_PROC_UPDATE_INTERVAL)
        except:
            if not been_warned[_WARN_IO_VAL]:
                been_warned[_WARN_IO_VAL] = True
                warn_q.put(_WARN_IO_VAL)
            continue
    # add last number of records reported from this process
    prog_q.put(num_recs_proc % _PROC_UPDATE_INTERVAL)

    return




def annotate_reads_with_fastq_main(args):
    fastq_slot = '/'.join(('/Analyses', args.basecall_group,
                           args.basecall_subgroup))
    if VERBOSE: th.status_message('Annotating FAST5s with sequence from FASTQs.')
    prog_q = Queue()
    warn_q = Queue()

    # open fast5 annotation processes
    ann_args = (single_fast5_q, fastq_recs, fastq_slot,
                prog_q, warn_q, bc_grp_name, bc_subgrp_name, overwrite)
    ann_ps = []
    for p_id in range(args.num_processes):
        ann_p = Process(target=_annotate_with_fastqs_worker, args=ann_args)
        ann_p.daemon = True
        ann_p.start()
        ann_ps.append(ann_p)

    main_wp_conn, wp_conn = Pipe()
    warn_prog_p = Process(target=_get_ann_queues,
                          args=(prog_q, warn_q, len(fast5s_read_ids), wp_conn))
    warn_prog_p.daemon = True
    warn_prog_p.start()

    fq_feed_p.join()
    for ann_p in ann_ps:
        ann_p.join()
    # send signal to warn/progress queue that all other processes are complete
    main_wp_conn.send(True)
    warn_prog_p.join()

    return













