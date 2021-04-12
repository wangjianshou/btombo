import sys
import io
import queue
import h5py
from multiprocessing import Process, Queue, Pipe
from time import sleep
from tqdm import tqdm
from tombo import tombo_helper as th
from tombo._preprocess import _prep_fast5_for_fastq # 获取single_fast5的read_id,创建注释slot 
from concurrent.futures import ThreadPoolExecutor, as_completed
from pudb.remote import set_trace

if sys.version_info[0] > 2:
        unicode = str

VERBOSE = False
#VERBOSE = True

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


def _get_ann_queues(prog_q, warn_q, wp_conn):
    #set_trace()
    if VERBOSE: bar = tqdm(smoothing=0)
    been_warned = dict((warn_code, False) for warn_code in _WARN_CODES)

    def update_warn(warn_val):
        if warn_val == _WARN_ID_VAL:
            if VERBOSE and not been_warned[_WARN_ID_VAL]:
                bar.write(
                    _WARN_PREFIX + 'Some FASTQ records contain read ' +
                    'identifiers not found in any FAST5 files or ' +
                    'sequencing summary files.',
                    file=sys.stderr)
            been_warned[_WARN_ID_VAL] = True
        elif warn_val == _WARN_IO_VAL:
            if VERBOSE and not been_warned[_WARN_IO_VAL]:
                bar.write(
                    _WARN_PREFIX + 'Some read files could not be accessed.',
                    file=sys.stderr)
            been_warned[_WARN_IO_VAL] = True
        elif warn_val == _WARN_MISMATCH_VAL:
            if VERBOSE and not been_warned[_WARN_MISMATCH_VAL]:
                bar.write(
                    _WARN_PREFIX + 'Read ID(s) found in sequencing summary ' +
                    'and FAST5 file are discordant. Skipping such reads.',
                    file=sys.stderr)
            been_warned[_WARN_MISMATCH_VAL] = True
        elif warn_val == _WARN_OVRWRT_VAL:
            if VERBOSE and not been_warned[_WARN_OVRWRT_VAL]:
                bar.write(
                    _WARN_PREFIX + 'Basecalls exsit in specified slot for ' +
                    'some reads. Set --overwrite option to overwrite these ' +
                    'basecalls.', file=sys.stderr)
            been_warned[_WARN_OVRWRT_VAL] = True
        else:
            if VERBOSE:
                bar.write('{}Invalid warning code encountered: {}'.format(
                    _WARN_PREFIX, warn_val), file=sys.stderr)

        return

    total_added_seqs = 0
    while True:
        try:
            iter_added = prog_q.get(block=False)
            total_added_seqs += iter_added
            if VERBOSE: bar.update(iter_added)
        except queue.Empty:
            try:
                warn_val = warn_q.get(block=False)
                update_warn(warn_val)
            except queue.Empty:
                sleep(0.1)
                # check if main thread has finished with all fastq records
                if wp_conn.poll():
                    break

    # collect all remaining warn and progress values
    while not prog_q.empty():
        iter_added = prog_q.get(block=False)
        total_added_seqs += iter_added
        if VERBOSE: bar.update(iter_added)
    while not warn_q.empty():
        warn_val = warn_q.get(block=False)
        update_warn(warn_val)
    if VERBOSE:
        bar.close()
        th.status_message('Added sequences to a total of ' +
                          str(total_added_seqs) + ' reads.')
        if total_added_seqs < num_read_ids:
            th.warning_message(
                'Not all read ids from FAST5s or sequencing summary files ' +
                'were found in FASTQs.\n\t\tThis can result from reads that ' +
                'failed basecalling or if full sets of FAST5s/sequence ' +
                'summaries are not processed with full sets of FASTQs.')

    return




def _annotate_with_fastqs_worker(single_fast5_q, anno_fast5_q, fastq_recs, fastq_slot,  
            prog_q, warn_q, bc_grp_name, bc_subgrp_name, overwrite): 
    been_warned = dict((warn_code, False) for warn_code in _WARN_CODES)
    num_recs_proc = 0
    while True:
        single_fast5 = single_fast5_q.get()
        #set_trace()
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
                    'Fastq', data=''.join(fastq_recs[file_parsed_id]),
                    dtype=h5py.special_dtype(vlen=unicode))
                num_recs_proc += 1
                if num_recs_proc % _PROC_UPDATE_INTERVAL == 0:
                    prog_q.put(_PROC_UPDATE_INTERVAL)
            #with io.open(single_fast5.name, 'wb') as f:
            #    f.write(single_fast5.getvalue())
            #single_fast5.close()
            anno_fast5_q.put(single_fast5)

        except:
            if not been_warned[_WARN_IO_VAL]:
                been_warned[_WARN_IO_VAL] = True
                warn_q.put(_WARN_IO_VAL)
            continue
    # add last number of records reported from this process
    prog_q.put(num_recs_proc % _PROC_UPDATE_INTERVAL)

    return




def annotate_reads_with_fastq_main(single_fast5_q, anno_fast5_q, fastq_fns, basecall_group,
                                   basecall_subgroup, overwrite, processes):
    fastq_slot = '/'.join(('/Analyses', basecall_group,
                           basecall_subgroup))
    if VERBOSE: th.status_message('Annotating FAST5s with sequence from FASTQs.')
    prog_q = Queue()
    warn_q = Queue()
    fastq_recs = {}
    for i in get_seq_recs_concurrent(fastq_fns):
        fastq_recs.update(i)
    # open fast5 annotation processes
    ann_args = (single_fast5_q, fastq_recs, fastq_slot,
                prog_q, warn_q, bc_grp_name, bc_subgrp_name, overwrite)
    ann_ps = []
    for p_id in range(processes):
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













