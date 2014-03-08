#!/usr/bin/env python

# standard library
import glob
import re
import shutil
from textwrap import dedent

# treeCl
from external import ExternalSoftware
from ..constants import TMPDIR
from ..datastructs.trcl_seq import TrClSeq
from ..errors import filecheck, directorymake, directoryquit, \
    optioncheck, OptionError, directorycheck
from ..utils.gapmasker import GapMasker
from ..utils import phymlIO


class Params(object):

    """ Class for writing parameter files for use with alfsim """

    def __init__(
        self,
        simulation_name='sim',
        working_directory='./alftmp/',
        outfile_path='./',
        unit_is_pam=True,
        ):

        if unit_is_pam is True:
            unit_is_pam = 'true'
        else:
            unit_is_pam = 'false'

        file_string = \
            dedent('''\
            # Filepath parameters
            ##############################
            # directories for file storage
            wdir := '{0}';
            dbdir := 'DB/';
            dbAncdir := 'DBancestral/';

            '''.format(working_directory))

        pam_string = \
            dedent('''\
            # Time units
            ##########################
            # timescale for simulation
            unitIsPam := {0}:

            '''.format(unit_is_pam))

        name_string = \
            dedent('''\
            # Simulation Name
            ####################
            # name of simulation
            mname := {0};

            '''.format(simulation_name))

        self.outfile_path = outfile_path
        self.name         = simulation_name

        self.files        = file_string
        self.pam          = pam_string
        self.name_string  = name_string
        self.genome       = ''
        self.subst        = ''
        self.indel_string = ''
        self.ratevar      = ''
        self.tree         = ''

    def __str__(self):

        sections = []
        for section in [
            'files',
            'pam',
            'name_string',
            'genome',
            'subst',
            'indel_string',
            'ratevar',
            'tree',
            ]:
            if hasattr(self, section):
                sections.append(getattr(self, section))
        return ''.join(sections)

    def filepaths(self, simulation_name='sim', working_directory='./'):

        file_string = \
            dedent('''\
            # Filepath parameters
            ##############################
            # directories for file storage
            wdir := '{0}';
            dbdir := 'DB/';
            dbAncdir := 'DBancestral/';

            # name of simulation
            mname := {1};

            '''.format(working_directory,
                   simulation_name))

        self.files = file_string
        return file_string

    def root_genome(
        self,
        number_of_genes=100,
        min_length=10,
        kappa=1,
        theta=1,
        ):

        genome_string = \
            dedent('''\
            # Root genome parameters
            #############################
            realseed := false;
            protStart := {0};
            minGeneLength := {1};
            gammaLengthDist := [{2},{3}];
            blocksize := 1:

            '''.format(number_of_genes,
                   min_length, kappa, theta))

        self.genome = genome_string
        return genome_string

    def pam(self, unit_is_pam=True):

        if unit_is_pam is True:
            unit_is_pam = 'true'
        else:
            unit_is_pam = 'false'

        pam_string = \
            dedent('''\
            # Time units
            ##########################
            # timescale for simulation
            unitIsPam := {0}:

            '''.format(unit_is_pam))

        self.pam = pam_string
        return pam_string

    def rename(self, name):

        name_string = \
            dedent('''\
            # Simulation Name
            ####################
            # name of simulation
            mname := {0};

            '''.format(name))

        self.name_string = name_string
        self.name = name
        return name_string

    def gtr_model(
        self,
        CtoT,
        AtoT,
        GtoT,
        AtoC,
        CtoG,
        AtoG,
        Afreq,
        Cfreq,
        Gfreq,
        Tfreq,
        allow_nonsense=False,
        ):

        if allow_nonsense:
            allow_nonsense = 'true'
        else:
            allow_nonsense = 'false'

        subst_string = \
            dedent('''\
            # Substitution Model Parameters
            ######################################################################################################
            substModels := [SubstitutionModel('GTR', [{0}, {1}, {2}, {3}, {4}, {5}], [{6}, {7}, {8}, {9}], {10})];

            '''.format(
            CtoT,
            AtoT,
            GtoT,
            AtoC,
            CtoG,
            AtoG,
            Afreq,
            Cfreq,
            Gfreq,
            Tfreq,
            allow_nonsense,
            ))

        self.subst = subst_string
        return subst_string

    def hky_model(
        self,
        alpha,
        beta,
        Afreq,
        Cfreq,
        Gfreq,
        Tfreq,
        allow_nonsense=False,
        ):

        if allow_nonsense:
            allow_nonsense = 'true'
        else:
            allow_nonsense = 'false'

        subst_string = \
            dedent('''\
            # Substitution Model Parameters
            #################################################################################
            substModels := [SubstitutionModel('HKY', [{0}, {1}], [{2}, {3}, {4}, {5}], {6})];

            '''.format(
            alpha,
            beta,
            Afreq,
            Cfreq,
            Gfreq,
            Tfreq,
            allow_nonsense,
            ))

        self.subst = subst_string
        return subst_string

    def f84_model(
        self,
        kappa,
        beta,
        Afreq,
        Cfreq,
        Gfreq,
        Tfreq,
        allow_nonsense=False,
        ):

        if allow_nonsense:
            allow_nonsense = 'true'
        else:
            allow_nonsense = 'false'

        subst_string = \
            dedent('''\
            # Substitution Model Parameters
            #################################################################################
            substModels := [SubstitutionModel('F84', [{0}, {1}], [{2}, {3}, {4}, {5}], {6})];

            '''.format(
            kappa,
            beta,
            Afreq,
            Cfreq,
            Gfreq,
            Tfreq,
            allow_nonsense,
            ))
        return subst_string

    def jc_model(self, allow_nonsense=False):

        if allow_nonsense:
            allow_nonsense = 'true'
        else:
            allow_nonsense = 'false'

        subst_string = \
            dedent('''\
            # Substitution Model Parameters
            ######################################################################
            substModels := [SubstitutionModel('F84', [1, 1], [seq(0.25,4)], {0})];

            '''.format(allow_nonsense))

        self.subst = subst_string
        return subst_string

    def one_word_model(self, model='WAG'):
        """ "One-word" models don't need any extra parameters. Codon: CPAM, ECM,
        ECMu AA:    GCB,  JTT, LG,  WAG """

        try:
            optioncheck(model, [
                'CPAM',
                'ECM',
                'ECMu',
                'GCB',
                'JTT',
                'LG',
                'WAG',
                ])
        except OptionError, e:
            print e
            return

        subst_string = \
            dedent('''\
            # Substitution Model Parameters
            ##########################################
            substModels := [SubstitutionModel('{0}')];

            '''.format(model))

        self.subst = subst_string
        return subst_string

    def indels(
        self,
        gain_rate=0.00002,
        gain_model='ZIPF',
        gain_params=[1.821],
        max_gain_length=10,
        loss_rate=0.00002,
        loss_model='ZIPF',
        loss_params=[1.821],
        max_loss_length=10,
        ):

        indel_string = \
            dedent('''\
            # Indel Parameters
            ####################################################################
            indelModels := [IndelModel({0}, {1}, {2}, {3}, {4}, {5}, {6}, {7})];

            '''.format(
            gain_rate,
            gain_model,
            gain_params,
            max_gain_length,
            loss_rate,
            loss_model,
            loss_params,
            max_loss_length,
            ))

        self.indel_string = indel_string
        return indel_string

    def rate_variation(
        self,
        shape=1,
        ncat=4,
        pinvar=0,
        ):
        """ Models rate variation among sites Only uses gamma distribution (ALF
        allows more options) """

        rate_var_string = \
            dedent('''\
            # Rate Variation Parameters
            ########################################################
            rateVarModels := [RateVarModel(Gamma, {0}, {1}, {2})];

            '''.format(ncat,
                   pinvar, shape))

        self.ratevar = rate_var_string
        return rate_var_string

    def custom_tree(self, treefile):

        custom_tree_string = \
            dedent('''\
            # Tree Parameters
            #####################
            treeType := 'Custom';
            treeFile := '{0}';

            '''.format(treefile))

        self.tree = custom_tree_string
        return custom_tree_string

    def BDtree(
        self,
        birthrate=0.01,
        deathrate=0.001,
        nspecies=15,
        mutation_rate=250,
        scale_tree=True,
        ultrametric=False,
        ):

        if ultrametric:
            ultrametric = 'true'
        else:
            ultrametric = 'false'

        if scale_tree:
            scale_tree = 'true'
        else:
            scale_tree = 'false'

        bdtree_string = \
            dedent('''\
            # Tree Parameters
            #####################
            treeType := 'BDTree';
            birthRate := {0};
            deathRate := {1};
            NSpecies := {2};
            ultrametric := {3};
            mutRate := {4};
            scaleTree := {5};

            '''.format(
            birthrate,
            deathrate,
            nspecies,
            ultrametric,
            mutation_rate,
            scale_tree,
            ))

        self.tree = bdtree_string
        return bdtree_string

    def write_parameters(self):
        outfile_name = '{0}/{1}.drw'.format(self.outfile_path, self.name)
        writer = open(outfile_name, 'w')
        writer.write(str(self))
        writer.flush()
        writer.close()

        return outfile_name


class ALF(ExternalSoftware):

    """ This class is a front end to the ALF simulator """

    default_binary = 'alfsim'

    def __init__(
        self,
        tree,
        datatype,
        tmpdir,
        num_genes=1,
        seqlength=10,
        gene_length_kappa=1,
        gene_length_theta=1,
        name='no_name',
        **kwargs
        ):

        super(ALF, self).__init__(tmpdir, **kwargs)
        self.tree = tree
        self.name = name
        self.num_genes = num_genes
        self.seqlength = seqlength
        self.datatype = optioncheck(datatype, ['protein', 'dna'])
        self.param_dir = \
            directorymake('{0}/alfsim_param_dir'.format(self.tmpdir))
        self.working_dir = \
            directorymake('{0}/alfsim_working_dir'.format(self.tmpdir))

        if datatype == 'dna':               # Length correction as ALF assumes
            seqlength = (seqlength / 3) + 1 # we mean amino acid sequences

        self.params = Params(simulation_name=name,
                             working_directory=self.working_dir,
                             outfile_path=self.param_dir)
        self.params.custom_tree('{0}/{1}.nwk'.format(self.param_dir, self.name))
        self.params.root_genome(number_of_genes=num_genes,
            min_length=seqlength, theta=gene_length_theta,
            kappa=gene_length_kappa)

    def __str__(self):
        return '\n'.join([
            '========================================',
            'ALF: alfsim sequence evolution simulator',
            '========================================\n',
            'Sim name: {0}'.format(self.name),
            'Temp dir: {0}'.format(self.tmpdir),
            'Datatype: {0}\n'.format(self.datatype),
            '----------------------------------------\n',
            str(self.params),
            ])

    def clean(self):
        shutil.rmtree(self.param_dir)
        shutil.rmtree(self.working_dir)

    def fix_alignment(
        self,
        alignment_file,
        replacement_dict,
        length=None,
        ):
        """ Alf uses its own names for taxa even when we provide a custom tree
        with our own taxon labels, so this function reapplies our names. At the
        same time we check that the alignments are the correct length """

        record = TrClSeq(alignment_file)
        record.headers = [replacement_dict[x[:x.rindex('/')]] for x in
                          record.headers]
        if length is not None:
            record.sequences = [seq[:length] for seq in record.sequences]

        record._update()
        record.sort_by_name()
        return record

    def read(self, verbosity=0, length_is_strict=False):
        alf_tree_file = \
            filecheck('{0}/{1}/RealTree.nwk'.format(self.working_dir,
                                self.name))
        alf_tree = open(alf_tree_file).read()
        replacement_dict = dict(zip(re.findall(r'(\w+)(?=:)', alf_tree),
                                re.findall(r'(\w+)(?=:)', self.tree.newick)))  # bug correction
        if verbosity > 0:
            print 'replacement_dict =', replacement_dict
        if self.datatype == 'dna':  # !!! ALF doesn't always write DNA
                                    # alignments
            alignments = \
                glob.glob('{0}/{1}/MSA/MSA_*_dna.fa'.format(self.working_dir,
                                    self.name))
        else:
            alignments = \
                glob.glob('{0}/{1}/MSA/MSA_*_aa.fa'.format(self.working_dir,
                                    self.name))
        if verbosity > 1:
            print alignments

        recs = []
        for f in alignments:
            if length_is_strict:
                rec = self.fix_alignment(f, replacement_dict, self.seqlength)
            else:
                rec = self.fix_alignment(f, replacement_dict)
            rec.datatype = self.datatype
            recs.append(rec)
        return recs

    def run(self, verbosity=0, cleanup=True, length_is_strict=False):
        params = self.write()
        if verbosity > 0:
            print 'Running ALF on {0}'.format(params)
        (stdout, stderr) = self.call()
        if verbosity > 1:
            print stdout, stderr
        records = self.read(verbosity, length_is_strict=length_is_strict)
        if cleanup:
            self.clean()
        return records

    def write(self):
        directorymake(self.param_dir)
        directorymake(self.working_dir)
        params = self.params.write_parameters()
        self.tree.write_to_file('{0}/{1}.nwk'.format(
            self.param_dir, self.name), internal_labels=False, scale=100)
        self.add_flag('', params)
        return params


# RUNNERS

def simulate_from_tree(
    tree,
    length,
    datatype='protein',
    tmpdir=TMPDIR,
    model='WAG',
    split_lengths=None,
    gene_names=None,
    ):

    directoryquit(tmpdir)
    optioncheck(datatype, ['dna', 'protein'])

    if datatype == 'dna':
        optioncheck(model, ['CPAM', 'ECM', 'ECMu'])
    else:
        optioncheck(model, [
                'CPAM',
                'ECM',
                'ECMu',
                'GCB',
                'JTT',
                'LG',
                'WAG',
                ])

    try:
        gamma = phymlIO.extract_gamma_parameter(tree)
    except:
        gamma = 1

    sim = ALF(tree, datatype, tmpdir, seqlength=length)

    sim.params.rate_variation(gamma)
    sim.params.one_word_model(model)

    record = sim.run()

    if split_lengths and gene_names:
        return record.split_by_lengths(split_lengths, gene_names)
    return record


def simulate_from_record(
    record,
    length=None,
    tmpdir=TMPDIR,
    model='WAG',
    allow_nonsense=True,
    split_lengths=None,
    gene_names=None,
    mask_gaps=True
    ):

    if not record.tree:
        raise Exception('Simulation template must have an associated tree')

    directorycheck(tmpdir)
    optioncheck(record.datatype, ['dna', 'protein'])

    if record.datatype == 'dna':
        optioncheck(model, ['CPAM', 'ECM', 'ECMu', 'GTR'])
    else:
        optioncheck(model, [
                'CPAM',
                'ECM',
                'ECMu',
                'GCB',
                'JTT',
                'LG',
                'WAG',
                ])

    length = (length if length is not None else record.seqlength)
    gamma = phymlIO.extract_gamma_parameter(record.tree)
    GTR_parameters = None

    sim = ALF(record.tree, record.datatype, tmpdir, seqlength=length)
    if model == 'GTR':
        GTR_parameters = phymlIO.extract_GTR_parameters(record.tree)
        sim.params.gtr_model(
        CtoT=GTR_parameters['CtoT'],
        AtoT=GTR_parameters['AtoT'],
        GtoT=GTR_parameters['GtoT'],
        AtoC=GTR_parameters['AtoC'],
        CtoG=GTR_parameters['CtoG'],
        AtoG=GTR_parameters['AtoG'],
        Afreq=GTR_parameters['Afreq'],
        Cfreq=GTR_parameters['Cfreq'],
        Gfreq=GTR_parameters['Gfreq'],
        Tfreq=GTR_parameters['Tfreq'],
        allow_nonsense=allow_nonsense,
        )
    else:
        sim.params.one_word_model(model)
    sim.params.rate_variation(gamma)

    sim_record = sim.run(length_is_strict=True)[0]
    print 'sim record: ', repr(sim_record)
    sim_record.tree = record.tree
    sim_record.name = record.name
    sim_record.datatype = record.datatype

    if mask_gaps:
        gp = GapMasker(record)
        gp.mask(sim_record)

    if split_lengths and gene_names:
        return sim_record.split_by_lengths(split_lengths, gene_names)

    return sim_record

def simulate_from_record_GTR(
    record,
    output_dir,
    name='tempsim',
    tmpdir=TMPDIR,
    allow_nonsense=True,
    split_lengths=None,
    gene_names=None,
    ):
    """ For use with trees estimated by Phyml -m GTR """

    length = record.seqlength
    tree = record.tree
    directorycheck(tmpdir)
    GTR_parameters = phymlIO.extract_GTR_parameters(tree)
    gamma = phymlIO.extract_gamma_parameter(tree)

    sim = ALF(tree, 'dna', tmpdir, length, name)

    sim.params.indels()
    sim.params.rate_variation(gamma)
    sim.params.gtr_model(
        CtoT=GTR_parameters['CtoT'],
        AtoT=GTR_parameters['AtoT'],
        GtoT=GTR_parameters['GtoT'],
        AtoC=GTR_parameters['AtoC'],
        CtoG=GTR_parameters['CtoG'],
        AtoG=GTR_parameters['AtoG'],
        Afreq=GTR_parameters['Afreq'],
        Cfreq=GTR_parameters['Cfreq'],
        Gfreq=GTR_parameters['Gfreq'],
        Tfreq=GTR_parameters['Tfreq'],
        allow_nonsense=allow_nonsense,
        )

    record = sim.run()

    if split_lengths and gene_names:
        with open('{0}/trees.txt'.format(output_dir), 'a') as trf:
            trf.write('{0}\t{1}\n'.format('-'.join(gene_names), tree.newick))
        for rec in record.split_by_lengths(split_lengths, gene_names):
            rec.write_phylip('{0}/{1}.phy'.format(output_dir, rec.name))
    else:
        with open('{0}/trees.txt'.format(output_dir), 'a') as trf:
            trf.write('{0}\t{1}\n'.format(record.name, tree.newick))
        record.write_phylip('{0}/{1}.phy'.format(output_dir, name))
