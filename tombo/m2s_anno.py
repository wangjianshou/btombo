from argparse import ArgumentParser
from multiprocessing import Process, Pool
import multiprocessing as mp
import threading
import logging
import os,sys
from ont_fast5_api import __version__
from ont_fast5_api.conversion_tools.conversion_utils import get_fast5_file_list, get_progress_bar
from ont_fast5_api.fast5_file import EmptyFast5, Fast5FileTypeError
from ont_fast5_api.fast5_interface import check_file_type, MULTI_READ
from tombo.multi_fast5 import MultiFast5File ###wjs
from tombo.bmk import F5BytesIO  ###wjs###
import io


VERBOSE = False

from tombo.anno_fastq import get_seq_recs_concurrent,get_seq_worker,_get_ann_queues,_annotate_with_fastqs_worker



logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
exc_info = False


def batch_convert_multi_files_to_single(input_path, output_folder, processes, threads_per_proc,
        fastq_fns, basecall_group, basecall_subgroup, overwrite, recursive, follow_symlinks):
    pool = Pool(processes)
    file_list = get_fast5_file_list(input_path, recursive, follow_symlinks=follow_symlinks)
    pbar = get_progress_bar(len(file_list))
    single_fast5_q = mp.Manager().Queue(20000)
    def update(result):
        input_file = result[0]
        with open(os.path.join(output_folder, "filename_mapping.txt"), 'a') as output_table:
            for filename in result[1]:
                output_table.write("{}\t{}\n".format(input_file, filename))
        pbar.update(pbar.currval + 1)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    results_array = []
    for batch_num, filename in enumerate(file_list):
        results_array.append(pool.apply_async(convert_multi_to_single,
                                              args=(filename, output_folder,
                                                    str(batch_num), single_fast5_q),
                                              callback=update))

    pool.close()


    fastq_slot = '/'.join(('/Analyses', basecall_group,
                           basecall_subgroup))
    if VERBOSE: th.status_message('Annotating FAST5s with sequence from FASTQs.')
    anno_fast5_q = mp.Queue(20000)
    prog_q = mp.Queue()
    warn_q = mp.Queue()
    fastq_recs = {}
    for i in get_seq_recs_concurrent(fastq_fns):
        fastq_recs.update(i)

    # open fast5 annotation processes
    ann_args = (single_fast5_q, anno_fast5_q, fastq_recs, fastq_slot, prog_q, warn_q,
                basecall_group, basecall_subgroup, overwrite)
    ann_ps = []

    #for p_id in range(threads_per_proc):
    for p_id in range(processes):
        ann_p = Process(target=_annotate_with_fastqs_worker, args=ann_args)
        #ann_p = threading.Thread(target=_annotate_with_fastqs_worker, args=ann_args)
        ann_p.daemon = True
        ann_p.start()
        ann_ps.append(ann_p)

    main_wp_conn, wp_conn = mp.Pipe()
    warn_prog_p = Process(target=_get_ann_queues,
                          args=(prog_q, warn_q, wp_conn))
    warn_prog_p.daemon = True
    warn_prog_p.start()

    write_threads = []
    for i in range(threads_per_proc):
        t = threading.Thread(target=single_fast5_write, args=(anno_fast5_q, ))
        t.start()
        write_threads.append(t)

    pool.join()

    for i in range(processes):
        single_fast5_q.put(None)
    for ann_p in ann_ps:
        ann_p.join()
    for i in range(threads_per_proc):
        anno_fast5_q.put(None)
    for t in write_threads:
        t.join()
    # send signal to warn/progress queue that all other processes are complete
    main_wp_conn.send(True)
    warn_prog_p.join()

    pbar.finish()


def single_fast5_write(single_fast5_q):
    while True:
        single_fast5 = single_fast5_q.get()
        if single_fast5 is None:
          break
        with io.open(single_fast5.name, 'wb') as f:
            f.write(single_fast5.getvalue())
        single_fast5.close()


def convert_multi_to_single(input_file, output_folder, subfolder, single_fast5_q):
    output_files = ()
    try:
        output_files = try_multi_to_single_conversion(input_file, output_folder, subfolder, single_fast5_q)
    except Exception as e:
        logger.error("{}\n\tFailed to copy files from: {}"
                     "".format(e, input_file), exc_info=exc_info)
    return input_file, output_files


def try_multi_to_single_conversion(input_file, output_folder, subfolder, single_fast5_q):
    output_files = []
    with MultiFast5File(input_file, 'r') as multi_f5:
        file_type = check_file_type(multi_f5)
        if file_type != MULTI_READ:
            raise Fast5FileTypeError("Could not convert Multi->Single for file type '{}' with path '{}'"
                                     "".format(file_type, input_file))
        for read in multi_f5.get_reads():
            try:
                output_file = os.path.join(output_folder, subfolder, "{}.fast5".format(read.read_id))
                create_single_f5(output_file, read, single_fast5_q)
                output_files.append(os.path.basename(output_file))
            except Exception as e:
                logger.error("{}\n\tFailed to copy read '{}' from {}"
                             "".format(str(e), read.read_id, input_file), exc_info=exc_info)
    return output_files


def create_single_f5(output_file, read, single_fast5_q):
    if not os.path.exists(os.path.dirname(output_file)):
        os.makedirs(os.path.dirname(output_file))
    output_file = F5BytesIO(output_file, True)  ###wjs###
    with EmptyFast5(output_file, 'w') as single_f5:
        for group in read.handle:
            if group == "Raw":
                read_number = read.handle["Raw"].attrs["read_number"]
                single_f5.handle.copy(read.handle[group], "Raw/Reads/Read_{}".format(read_number))
            elif group in ("channel_id", "context_tags", "tracking_id"):
                if "UniqueGlobalKey" not in single_f5.handle:
                    single_f5.handle.create_group("UniqueGlobalKey")
                single_f5.handle.copy(read.handle[group], "UniqueGlobalKey/{}".format(group))
            else:
                single_f5.handle.copy(read.handle[group], group)
    ###wjs###
    single_fast5_q.put(output_file)
    #with io.open(output_file.name, 'wb') as f:
    #    f.write(output_file.getvalue())
    #output_file.close()
    ###wjs###


def main():
    parser = ArgumentParser("")
    parser.add_argument('-i', '--input_path', required=True,
                        help="MultiRead fast5 file or path to directory of MultiRead files")
    parser.add_argument('-s', '--save_path', required=True,
                        help="Folder to output SingleRead fast5 files to")
    parser.add_argument('-fq', '--fastq-filenames', required=True, nargs='+',
                        help='FASTQ filenames containing basecalls to be added to ' +
                        'raw FAST5 files.')
    parser.add_argument('--recursive', action='store_true',
                        help="Search recursively through folders for MultiRead fast5 files")
    parser.add_argument('--ignore_symlinks', action='store_true',
                        help="Ignore symlinks when searching recursively for fast5 files")
    parser.add_argument('-p', '--process', type=int, default=1, required=False,
                        help="Number of threads to use")
    parser.add_argument('-t', '--threads_per_proc', type=int, default=6, required=False,
                        help="Number of threads to use")
    parser.add_argument('--basecall-group', default='Basecall_1D_000', 
                        help="FAST5 group obtain original basecalls (under Analyses group)."+
                             'Default: %(default)s')
    parser.add_argument('--basecall-subgroup', default='BaseCalled_template', 
                        help="FAST5 subgroup (under /Analyses/[--basecall-group]/) under which" +
                             'to store basecalls from FASTQs. Default: %(default)s')
    parser.add_argument('--overwrite', default=False, action='store_true', 
                        help='Overwrite previous corrected group in FAST5 files. Note: ' +
                        'only effects --corrected-group or --new-corrected-group.')
    

    parser.add_argument('-v', '--version', action='version', version=__version__)
    args = parser.parse_args()

    batch_convert_multi_files_to_single(args.input_path, args.save_path, args.process, args.threads_per_proc,
                                        args.fastq_filenames, args.basecall_group, args.basecall_subgroup,
                                        args.overwrite, args.recursive, follow_symlinks=not args.ignore_symlinks)


if __name__ == '__main__':
    main()
