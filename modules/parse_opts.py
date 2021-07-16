#!/usr/bin/env python
#title       : parse_options.py
#description : Parse the input parameters of the deconPeaker algorithm
#author      : Huamei Li
#date        : 03/06/2018
#type        : module
#version     : 2.7

#--------------------------------------------------------
# load own modules

from modules.utils    import *
from modules.opt_cmds import opts

#--------------------------------------------------------
# global setting

LOGS = log_infos() # logging informative

#--------------------------------------------------------

def read_snyaml(yaml_fil, ref=True):
    '''
    parse YAML file, and then convert it into pd.DataFrame data structure
    :param yaml_fil: [str] the path of YAML format file
    :param ref: [bool] read YAML file of pure samples or not
    :reuturn: sn_infos [pd.dataframe]
    
    '''
    keys = ['bam', 'peak', 'cellname', 'batch'] if ref else ['bam', 'label']
    sn_dict, sn_infos = __import__('yaml').load(open(yaml_fil)), []
    
    for sn, infos in sn_dict.items():
        infos = { key : value if isinstance(value, list) else [value] for key, value in infos.items() }
        if ref and 'batch' not in infos:
            tmp_name = itemgetter(*np.random.choice(51, 10))(__import__('string').letters)
            infos['batch'] = ['_XX_{}'.format(''.join(tmp_name))]
        try:
            sub_info = [value for key in keys for value in infos[key]]
        except KeyError:
            LOGS.error('YAML key error, {} only contains the following keys: {}'.format(yaml_fil, ', '.join(keys)))
            raise
        sn_infos.append(sub_info)

    col_names = ['BAM', 'PEAK', 'CELL', 'REP_STATUS'] if ref else ['BAM', 'CELL']
    sn_infos = pd.DataFrame(sn_infos, columns=col_names)
    return sn_infos

def die(func, msg):
    '''
    print help infos and die the program
    :param func: [function] excuted main program
    :param msg: [str] error message
    :return: 0
    
    '''
    func.print_help()
    sys.exit(LOGS.error(msg))
    return 0

def basic_parameters(func, cmds, values, file_pos):
    '''
    check whether basic parameters correct or not
    :param func: [function] excuted main program
    :param cmds: [list] command string
    :param values: [list] input values by users
    :param file_pos: [list] which value represents file
    :return: 0
    
    '''
    infos = zip(cmds, values)
    for idx, values in enumerate(infos):
        cmd, vv = values
        if not vv:
            die(func, 'Too few input parameters, {} must has input, exiting......'.format(cmd))
        if idx not in file_pos: continue
        if not os.path.exists(vv):
            die(func, 'File cannot be detected of {} option'.format(cmd))
    return 0

def yaml_keys(func, keys, yamlfile):
    '''
    check keys in YAML file or not
    :param func: [function] excuted main program
    :param keys: [list] keys which need to be check
    :param yamlfile: [str/file] YAML file
    :return: infos [dict] YAML infos
    
    '''
    infos = __import__('yaml').load(open(yamlfile))
    for key, value in infos.iteritems():
        for subkey in keys:
            if subkey in value:
                try:
                    boolv = isinstance(value[subkey][0], list)
                except TypeError:
                    continue
                if not boolv: continue
                infos[key][subkey] = infos[key][subkey][0]
                continue
            msg = 'YAML key error, {} only contains the following keys: {}, exiting......'.format(yamlfile, ', '.join(keys))
            die(func, msg)
    return infos

def preprocess():
    '''
    check whether the parameters of preprocess step are legal
    :return: 0
    
    '''
    if ARGS.pure: ARGS.infos = read_snyaml(ARGS.pure, ref=True)
    ARGS.blacklist = 'extra/{}.blacklist.bed'.format(ARGS.hg_version) #determine the genome version of blacklist
    ARGS.hg_genome = 'extra/{}_tss.pos'.format(ARGS.hg_version)

    ARGS.qvalue = -np.log10(ARGS.qvalue) # log transfer for q-value (log10Pvalue)
    if not ARGS.prefix:
        die(ARGS.preprocess, '--prefix must be assigned, exiting......')
    return 0

def findctsps():
    '''
    check whether the parameters of findctsps are legal
    :return: 0
    
    '''
    cmds, values = ['--profile/-f', '--phenotype/-c'], \
            [ARGS.profile, ARGS.phenotype]
    basic_parameters(ARGS.findctsps, cmds, values, [0, 1])
    ARGS.merge_replicates = 'mean'
    return 0

def deconvolution():
    '''
    check whether the parameters of deconvolution step are legal
    :return: 0
    
    '''
    cmds, values = ['--pure/-p', '--mixture/-m'], [ARGS.pure, ARGS.mixture]
    basic_parameters(ARGS.deconv, cmds, values, [0, 1])
    
    if not (ARGS.mixture and ARGS.format):
        die(ARGS.deconv, '--format/-f name must be assigned when --mixture/-m exist, exiting......')
    if ARGS.format == 'BAM':
        ARGS.infos, ARGS.mapq = read_snyaml(ARGS.mixture, ref=False), 10
    
    ARGS.pvalue = True if ARGS.pvalue == 'TRUE' else False
    return 0

def simulate():
    '''
    check whether the parameters of simulate step are legal
    :return: 0
    
    '''
    cmds, values, purekeys = ['--pure/-p'], [ARGS.pure], ['data', 'cellname']
    basic_parameters(ARGS.simulate, cmds, values, [0])
    
    if not ARGS.format:
        die(ARGS.simulate, '--format/-f must be assigned, optional value is BAM or BED, exiting......')
    if ARGS.format == 'BED' and not ARGS.genome:
        die(ARGS.simulate, '--genome/-g must be assigned when --format/-f=BED, exiting......')
    
    ARGS.pure_infos = dict2pd(yaml_keys(ARGS.simulate, purekeys, ARGS.pure), purekeys)
    pure_cellnum = len(ARGS.pure_infos.cellname.unique())
    ARGS.mixture = [ map(int, ARGS.mixture.strip().split(',')) if ARGS.mixture else range(1, pure_cellnum + 1) ][0]
    bool_v = np.sum(np.array(ARGS.mixture) <= 0)
    
    if bool_v:
        die(ARGS.simulate, '--mixture/-m assigned each value must be greater than 0, exiting......')
    return 0

def check_platform():
    '''
    check unix platform or not
    :return: 0

    '''
    platname = __import__('platform').system()
    if platname != 'Linux':
        sys.exit(LOGS.error('Does not support {}, only supports Linux series platform'.format(platname)))
    return 0

def parse_opts():
    '''
    parse all input parameters
    :return: ARGS [object] global variable which cointains total input parameters of deconPeaker
    
    '''
    global ARGS
    ARGS = opts()
    funcs, values = [preprocess, findctsps, deconvolution, simulate], \
                    ['preprocess', 'findctsps', 'deconvolution', 'simulation']
    
    if funcs == preprocess: check_platform()
    funcs[values.index(ARGS.sub_parser)]()
    
    outdir = os.path.join(ARGS.outdir, ARGS.sub_parser)
    status, ARGS.outdir = create_dirs(outdir), outdir
    [ARGS.__delattr__(attr) for attr in ['preprocess', 'findctsps', 'deconv', 'simulate']]
    return ARGS
