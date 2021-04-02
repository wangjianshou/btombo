from __future__ import division, unicode_literals, absolute_import

from builtins import int, range, dict, map, zip

import os
import io
import re
import sys
import queue
import traceback
import threading

import queue
from pudb.remote import set_trace

# pip allows tombo install without correct version of mappy, so check here
try:
    import mappy
    if sys.version_info[0] > 2:
        mappy.Aligner(os.path.devnull).seq('')
    else:
        mappy.Aligner(os.path.devnull).seq(b'')
except AttributeError:
    th.error_message_and_exit('Tombo requires mappy version >= 2.10.')

# Future warning from cython in h5py
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import h5py

import numpy as np
np.seterr(all='raise')
import multiprocessing as mp

from tqdm import tqdm
from tqdm._utils import _term_move_up

from time import sleep
from operator import itemgetter
from collections import defaultdict
from pkg_resources import resource_string


from . import bmk


if sys.version_info[0] > 2:
    unicode = str

# import tombo modules/functions
from . import tombo_stats as ts
from . import tombo_helper as th

from ._default_parameters import (
    EXTRA_SIG_FACTOR, MASK_FILL_Z_SCORE,
    MASK_BASES, DEL_FIX_WINDOW, MAX_DEL_FIX_WINDOW,
    MIN_EVENT_TO_SEQ_RATIO, MAX_RAW_CPTS, SHIFT_CHANGE_THRESH,
    SCALE_CHANGE_THRESH, SIG_MATCH_THRESH, DNA_SAMP_TYPE, RNA_SAMP_TYPE,
    USE_RNA_EVENT_SCALE, RNA_SCALE_NUM_EVENTS, RNA_SCALE_MAX_FRAC_EVENTS,
    START_CLIP_PARAMS, STALL_PARAMS, COLLAPSE_RNA_STALLS, COLLAPSE_DNA_STALLS)
START_CLIP_PARAMS = th.startClipParams(*START_CLIP_PARAMS)
DEFAULT_STALL_PARAMS = th.stallParams(**STALL_PARAMS)

from ._c_dynamic_programming import (
    c_reg_z_scores, c_banded_forward_pass, c_base_z_scores,
    c_base_forward_pass, c_base_traceback)


# list of classes/functions to include in API
__all__ = [
    'get_read_seq', 'map_read', 'resquiggle_read',
    'segment_signal', 'find_adaptive_base_assignment',
    'resolve_skipped_bases_with_raw', 'find_seq_start_in_events',
    'find_static_base_assignment']


VERBOSE = True

_PROFILE_RSQGL = False
_PROFILE_IO_MAP = False

# use (mapping) clipped bases at the start of read to identify start position
USE_START_CLIP_BASES = False

# experimental RNA adapter trimming
TRIM_RNA_ADAPTER = False


# text output debugging
_DEBUG_PARAMS = False
_DEBUG_BANDWIDTH = False
_DEBUG_START_BANDWIDTH = False

# plot output debugging
_DEBUG_DP_ENDS = False
_DEBUG_DP_START = False
_DEBUG_CLIP_START = False
# fit debug plot requires r cowplot package to be installed
_DEBUG_FIT = False
_DEBUG_START_CLIP_FIT = False
# raw signal re-squiggle DP
_DEBUG_RAW_DP = False

# don't plot more than one debug type at a time
assert sum((
    _DEBUG_DP_ENDS, _DEBUG_FIT, _DEBUG_START_CLIP_FIT,
    _DEBUG_DP_START, _DEBUG_CLIP_START, _DEBUG_RAW_DP)) <= 1
_DEBUG_PLOTTING = any((
    _DEBUG_FIT, _DEBUG_START_CLIP_FIT, _DEBUG_DP_ENDS, _DEBUG_DP_START,
    _DEBUG_CLIP_START, _DEBUG_RAW_DP))
_DRY_RUN = any((
    _DEBUG_PARAMS, _DEBUG_BANDWIDTH, _DEBUG_START_BANDWIDTH, _DEBUG_PLOTTING))

_UNEXPECTED_ERROR_FN = 'unexpected_tombo_errors.{}.err'
_MAX_NUM_UNEXP_ERRORS = 50
_MAX_QUEUE_SIZE = 1000

_BROKEN_PIPE_ERR = (
    'Connection to re-squiggle process broken. THREAD CANNOT BE RECOVERED.')


##################################
########## Debug Output ##########
##################################

from .resquiggle import _write_params_debug

from .resquiggle import _debug_plot_dp
from .resquiggle import _debug_raw_dp

from .resquiggle import _debug_fit
from .resquiggle import _open_debug_pdf
from .resquiggle import _close_debug_pdf


############################################
########## Raw Signal Re-squiggle ##########
############################################

from .resquiggle import raw_forward_pass
from .resquiggle import raw_traceback
from .resquiggle import resolve_skipped_bases_with_raw

#####################################################
########## Static Band Dynamic Programming ##########
#####################################################

from .resquiggle import find_static_base_assignment

#######################################################
########## Adaptive Band Dynamic Programming ##########
#######################################################

from .resquiggle import _get_masked_start_fwd_pass

from .resquiggle import find_seq_start_in_events

from .resquiggle import find_seq_start_from_clip_basecalls

from .resquiggle import get_rel_raw_coords

from .resquiggle import find_adaptive_base_assignment

######################################
########## Re-squiggle Read ##########
######################################

from .resquiggle import segment_signal

from .resquiggle import resquiggle_read


#######################################
########## Genomic Alignment ##########
#######################################

from .resquiggle import get_read_seq

from .resquiggle import map_read

from .resquiggle import _io_and_map_read

#########################################
########## Re-squiggle Workers ##########
#########################################

from .resquiggle import _resquiggle_worker

if _PROFILE_RSQGL:
    _resquiggle_wrapper = _resquiggle_worker
    def _resquiggle_worker(*args):
        import cProfile
        cProfile.runctx('_resquiggle_wrapper(*args)', globals(), locals(),
                        filename='resquiggle_main.prof')
        return

from .resquiggle import _io_and_mappy_thread_worker

if _PROFILE_IO_MAP:
    _io_map_wrapper = _io_and_mappy_thread_worker
    def _io_and_mappy_thread_worker(*args):
        import cProfile
        cProfile.runctx('_io_map_wrapper(*args)', globals(), locals(),
                        filename='io_map_main.prof')
        return


############################################
########## Multi-process Handling ##########
############################################

from .resquiggle import _get_progress_fail_queues

from .resquiggle import _get_index_queue

from .resquiggle import _fill_files_queue

#from .resquiggle import resquiggle_all_reads


###################################
########## Main Function ##########
###################################

from .resquiggle import _parse_files_and_lock_dirs

#from .resquiggle import _resquiggle_main

def resquiggle_all_reads(
        single_fast5_q, num_reads, aligner, bc_grp, bc_subgrps, corr_grp, std_ref,
        seq_samp_type, outlier_thresh, overwrite, num_ps, threads_per_proc,
        compute_sd, skip_index, rsqgl_params, save_params, sig_match_thresh,
        obs_filter, const_scale, q_score_thresh, skip_seq_scaling,
        max_scaling_iters, failed_reads_fn, fast5s_basedir, num_update_errors,
        sig_len_rng, seq_len_rng):
    """Perform genomic alignment and re-squiggle algorithm
    """
    failed_reads_q = mp.Queue()
    index_q = mp.Queue(maxsize=_MAX_QUEUE_SIZE) if not skip_index else None 
    progress_q = mp.Queue()

    # open all multiprocessing pipes and queues before threading
    # as opening threads before all process are open seems to cause
    # a deadlock when some processes are started.
    # starting all multiprocess objects seems to fix this.

    map_conns = [] 
    rsqgl_ps = [] 
    for _ in range(num_ps):
        proc_rsqgl_conns = [] 
        for _ in range(threads_per_proc):
            # open mp pipe to communicate with re-squiggle process
            map_conn, rsqgl_conn = mp.Pipe()
            map_conns.append(map_conn)
            proc_rsqgl_conns.append(rsqgl_conn)
        # open re-squiggle process to void intensive processing hitting the GIL
        rsqgl_args = (
            proc_rsqgl_conns, std_ref, outlier_thresh, corr_grp, seq_samp_type,
            rsqgl_params, save_params, const_scale, skip_seq_scaling,
            max_scaling_iters)
        rsqgl_process = mp.Process(target=_resquiggle_worker, args=rsqgl_args)
        rsqgl_process.daemon = True
        rsqgl_process.start()
        rsqgl_ps.append(rsqgl_process)

    # failed read and progress queues getter
    main_pf_conn, pf_conn = mp.Pipe()
    pf_p = mp.Process(target=_get_progress_fail_queues,
                      args=(progress_q, failed_reads_q, pf_conn, num_reads,
                            failed_reads_fn, num_update_errors))
    pf_p.daemon = True
    pf_p.start()

    # index queue getter
    if index_q is not None:
        main_index_conn, index_conn = mp.Pipe()
        index_p = mp.Process(target=_get_index_queue, args=(
            index_q, index_conn, fast5s_basedir, corr_grp))
        index_p.daemon = True
        index_p.start()
    # now open mapping thread for each map connection created above
    map_and_io_ts = []
    for map_conn in map_conns:
        map_args = (single_fast5_q, progress_q, failed_reads_q, index_q, bc_grp,
                    bc_subgrps, corr_grp, aligner, outlier_thresh, compute_sd,
                    sig_match_thresh, obs_filter, seq_samp_type, overwrite,
                    map_conn, q_score_thresh, std_ref, sig_len_rng, seq_len_rng)
        t = threading.Thread(target=_io_and_mappy_thread_worker,
                             args=map_args)
        t.daemon = True
        t.start()
        map_and_io_ts.append(t)

    # wait for all mapping and re-squiggling workers to finish
    for t in map_and_io_ts:
        t.join()
    for rsqgl_p in rsqgl_ps:
        rsqgl_p.join()

    # in a very unlikely case the progress/fail queue could die while the
    # main process remains active and thus we would have a deadlock here
    if pf_p.is_alive():
        # send signal to getter queue to finish and return results
        main_pf_conn.send(True)
        # returns total number of processed reads if that is needed
        main_pf_conn.recv()
        pf_p.join()
    if index_q is not None:
        main_index_conn.send(True)
        index_p.join()

    return





def pipeline_main(args):
    """Main method for pipeline
    """
    if args.processes > 1 and _DEBUG_PLOTTING:
        th.error_message_and_exit(
            'Cannot run multiple processes and debugging.')
    if _DEBUG_PLOTTING:
        th.warning_message(
            'Producing de-bug plotting output. Can be very slow and should ' +
            'only be run on a small number of files.')
    if _DRY_RUN:
        th.warning_message(
            'Producing de-bug output. Not saving re-squiggle results.')
    if _DEBUG_BANDWIDTH or _DEBUG_START_BANDWIDTH:
        sys.stdout.write(
            'bandwidth\tmin_bw_edge_buffer\tmean_dp_score\tread_id\n')

    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE
    ts.VERBOSE = VERBOSE

    if args.print_advanced_arguments:
        from . import _option_parsers
        _option_parsers.print_advanced_resquiggle()
        sys.exit()

    if args.basecall_group == args.corrected_group:
        th.error_message_and_exit(
            '--basecall-group and --corrected-group must ' +
            'be different.')
    # check simple arguments for validity first
    outlier_thresh = args.outlier_threshold
    if outlier_thresh is not None and outlier_thresh <= 0:
        outlier_thresh = None
    obs_filter = th.parse_obs_filter(args.obs_per_base_filter) \
                 if 'obs_per_base_filter' in args else None

    if VERBOSE: th.status_message('Loading minimap2 reference.')
    # to be enabled when mappy genome sequence extraction bug is fixed
    aligner = mappy.Aligner(
        str(args.reference), preset=str('map-ont'), best_n=1)
    if not aligner:
        th.error_message_and_exit(
            'Failed to load reference genome FASTA for mapping.')

    # get files as late as possible in startup since it takes the longest
    # and so other errors can't happen after locks are written
    files, fast5s_basedir, lock_fns = _parse_files_and_lock_dirs(args)


    ###wjs###
    num_ps = args.processes
    threads_per_proc = args.threads_per_process
    single_fast5_q = queue.Queue(_MAX_QUEUE_SIZE)
    single_fast5 = (bmk.F5BytesIO(i) for i in files)
    num_reads = len(files)
    single_fast5_t = threading.Thread(target=_fill_files_queue,
                                      args=(single_fast5_q, single_fast5, num_ps * threads_per_proc))
    single_fast5_t.start()

    ###wjs###

    try:
        seq_samp_type = None
        if args.seq_sample_type is not None:
            seq_samp_type = th.seqSampleType(RNA_SAMP_TYPE, True) \
                            if args.seq_sample_type == RNA_SAMP_TYPE else \
                               th.seqSampleType(DNA_SAMP_TYPE, False)
        if args.tombo_model_filename is not None and seq_samp_type is None:
            seq_samp_type = th.get_seq_sample_type(fast5_fns=files)
        # parse tombo model
        std_ref = ts.TomboModel(
            ref_fn=args.tombo_model_filename, seq_samp_type=seq_samp_type,
            fast5_fns=files)
        seq_samp_type = std_ref.seq_samp_type
        if seq_samp_type.name == DNA_SAMP_TYPE and USE_START_CLIP_BASES:
            std_ref = std_ref.reverse_sequence_copy()

        sig_match_thresh = args.signal_matching_score
        if sig_match_thresh is None:
            sig_match_thresh = SIG_MATCH_THRESH[seq_samp_type.name]

        const_scale = None
        if args.fixed_scale is not None:
            const_scale = args.fixed_scale
        elif args.fit_global_scale:
            const_scale = ts.estimate_global_scale(files)

        rsqgl_params = ts.load_resquiggle_parameters(
            seq_samp_type, args.signal_align_parameters,
            args.segmentation_parameters)
        save_params = ts.load_resquiggle_parameters(
            seq_samp_type, args.signal_align_parameters,
            args.segmentation_parameters, use_save_bandwidth=True)


        resquiggle_all_reads(
            single_fast5_q, num_reads, aligner, args.basecall_group, args.basecall_subgroups,
            args.corrected_group, std_ref, seq_samp_type, outlier_thresh,
            args.overwrite, args.processes, args.threads_per_process,
            args.include_event_stdev, args.skip_index,
            rsqgl_params, save_params, sig_match_thresh,
            obs_filter, const_scale, args.q_score,
            args.skip_sequence_rescaling, args.max_scaling_iterations,
            args.failed_reads_filename, fast5s_basedir,
            args.num_most_common_errors, args.signal_length_range,
            args.sequence_length_range)

    finally:
        th.clear_tombo_locks(lock_fns)
        #single_fast5_t.join()
    return




if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `tombo -h`')
    sys.exit(1)
