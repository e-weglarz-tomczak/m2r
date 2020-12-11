import os
import pickle
import cobra

class Microbiota2Recon(object):
    def __init__(self):
        print('\n\nWelcome to Microbiota-to-Recon (M2R) program!\n')
        print(' ***************************\n',
              '* __   __   ___     ____  *\n',
              '*|  \ /  | / _ \   |  _ \ *\n',
              '*|   ^   | \/ \ |  | |_|/ *\n',
              '*|  | |  |   / /_  |   \  *\n',
              '*|__| |__|  /____\ |_|\_\ *\n',
              '*                         *\n',
              '***************************')
        print('© E. Weglarz-Tomczak & J. Tomczak\n\nThis program consists of three steps that modify RECON\n'
              'by introducing information from microbiota.\n')

        print('------------------\n')

        self.folder = input('Please provide a folder name where the microbiota files (.mat) are (default: microbiota): ') or 'microbiota'
        print('\t -> {}'.format(self.folder))

        self.files = [item for item in input('Please provide a list of names of microbiota that are supposed to influence the model separated by a space (e.g., Ruminococcus_sp_5_1_39BFAA.mat Acinetobacter_baumannii_AB0057.mat, default: load all files from the folder): ').split()]
        if len(self.files) == 0:
            print('\t -> all files in the directory will be loaded')
        else:
            print('\t -> {}'.format(self.files))

        while True:
            self.alpha = float(input('Please provide a value of alpha (between 0 and 1; default: 0.5): ') or "0.5")
            if self.alpha <= 0. or self.alpha >= 1.:
                print("Sorry, alpha must be between 0. and 1. Please provide a new value.")
                continue
            else:
                break
        print('\t -> {}'.format(self.alpha))

        while True:
            self.fin_value = float(input('Please provide a value that fluxes should be normalized to (default: 1000.): ') or "1000")
            if self.fin_value <= 0.:
                print("Sorry, the normalization values must larger than 0. Please provide a new value.")
                continue
            else:
                break
        print('\t -> {}'.format(self.fin_value))

        self.recon_dir = input('Please provide a directory name where the model is (default: models): ') or 'models'
        print('\t -> {}'.format(self.recon_dir))

        self.recon_file = input('Please provide a name of the model (please rememeber about the extension: either .mat, .json, .xml, .yaml or .pickle): ')
        print('\t -> {}'.format(self.recon_file))

        print('\t Loading the model (this may take a while)...')
        if self.recon_file[-3:] == 'mat':
            self.R = cobra.io.load_matlab_model(os.path.join(self.recon_dir, self.recon_file))
        elif self.recon_file[-3:] == 'son':
            self.R = cobra.io.load_json_model(os.path.join(self.recon_dir, self.recon_file))
        elif self.recon_file[-3:] == 'aml':
            self.R = cobra.io.load_yaml_model(os.path.join(self.recon_dir, self.recon_file))
        elif self.recon_file[-3:] == 'xml':
            self.R = cobra.io.read_sbml_model(os.path.join(self.recon_dir, self.recon_file))
        else:
            with open(os.path.join(self.recon_dir, self.recon_file), 'rb') as handle:
                self.R = pickle.load(handle)
        print('\t\t ...success!')

        self.names = []
        for i in range(len(self.R.metabolites)):
            metab = self.R.metabolites[i]
            self.names.append(metab.id)

        self.new_recon_name = input('Please provide a name for the new model (possible extensions: .mat, .json, .xml, .yaml or .pickle, e.g., New_recon.mat): ')
        print('\t -> {}'.format(self.new_recon_name))

        print('\nAll is set up!\nLet`s start!\n')
        print('------------------\n')

    @staticmethod
    def get_ins_outs(pf, inputs=None, outputs=None):
        if inputs is None:
            inputs = {}
        if outputs is None:
            outputs = {}

        for i in range(len(pf)):
            name = pf.iloc[i].name
            if pf.loc[name, 'is_input'] == 1.:
                if name in inputs.keys():
                    inputs[name] = inputs[name] + pf.loc[name, 'flux']
                else:
                    inputs[name] = pf.loc[name, 'flux']
            elif pf.loc[name, 'is_input'] == 0.:
                if name in outputs.keys():
                    outputs[name] = outputs[name] + pf.loc[name, 'flux']
                else:
                    outputs[name] = pf.loc[name, 'flux']
        return inputs, outputs

    def calculate_ins_outs(self, files, folder):
        inputs = None
        outputs = None

        # if files is empty, load all microbiota!
        if len(files) == 0:
            files = os.listdir(folder)

        print('Step 1/3: Read microbiota...')
        for i in range(len(files)):
            microb = files[i]
            print('\t(File {} out of {}) {}'.format(i+1, len(files), microb))
            # load microbiota
            if microb[-3:] == 'mat':
                M = cobra.io.load_matlab_model(os.path.join(folder, microb))
            else:
                with open(os.path.join(folder, microb), 'rb') as handle:
                    M = pickle.load(handle)
            # get summary to read fluxes and inputs and outputs
            pf = M.summary()._generate()
            # get inputs and outputs from the microbiota
            inputs, outputs = self.get_ins_outs(pf, inputs, outputs)

        print('\t...done!\n')

        return inputs, outputs

    def normalize_ins_outs(self, inputs, outputs, fin_value=1000.):
        print('Step 2/3:  Normalize metabolites...')
        max_val = max(max(list(inputs.values())), max(list(outputs.values())))
        # inputs
        for ki in inputs:
            inputs[ki] = fin_value * float(inputs[ki] / max_val)

        # outputs
        for ko in outputs:
            outputs[ko] = fin_value * float(outputs[ko] / max_val)

        print('\t... done!\n')
        return inputs, outputs

    def modify_recon(self, Recon_org, inputs, outputs, alpha=0.5):
        # alpha: how much microbiota is weighted in relation to RECON
        assert (alpha > 0.) and (alpha < 1.), "Weight must be between 0. and 1.!"

        Recon = Recon_org.copy()

        print('Step 3/3:  Modify the model...')
        # scale all lower bounds of the EXTERNAL reactions
        for r in Recon.reactions:
            if 'EX_' in r.id:
                r.lower_bound = alpha * r.lower_bound

        # go over inputs and modify lower bounds
        for ki in inputs:
            if 'EX_' + ki in self.names:
                Recon.reactions.get_by_id('EX_' + ki).lower_bound = Recon.reactions.get_by_id('EX_' + ki).lower_bound + (1. - alpha) * inputs[ki]

        # go over outputs
        for ko in outputs:
            if 'EX_' + ko in self.names:
                Recon.reactions.get_by_id('EX_' + ko).lower_bound = Recon.reactions.get_by_id('EX_' + ko).lower_bound - (1. - alpha) * outputs[ko]

        print('\t... done!\n')
        return Recon

    def procedure(self):
        ins, outs = self.calculate_ins_outs(self.files, self.folder)

        ins_n, outs_n = self.normalize_ins_outs(ins, outs)

        Recon = self.modify_recon(self.R, ins_n, outs_n, alpha=self.alpha)

        if self.new_recon_name[-3:] == 'mat':
            cobra.io.save_matlab_model(Recon, os.path.join(self.recon_dir, self.new_recon_name))
        elif self.new_recon_name[-3:] == 'son':
            cobra.io.save_json_model(Recon, os.path.join(self.recon_dir, self.new_recon_name))
        elif self.new_recon_name[-3:] == 'aml':
            cobra.io.save_yaml_model(Recon, os.path.join(self.recon_dir, self.new_recon_name))
        elif self.new_recon_name[-3:] == 'xml':
            cobra.io.write_sbml_model(Recon, os.path.join(self.recon_dir, self.new_recon_name))
        else:
            with open(os.path.join(self.recon_dir, self.new_recon_name), 'wb') as handle:
                pickle.dump(Recon, handle)

        print('\nThank you for using M2R!\n')
        print('------------------')
        print('   _____    _   _    ____   _\n',
              ' |  _  \  \ \ / /  |  __| | |\n',
              ' | |_| /   \   /   | |__  | |\n',
              ' | |_| \    \ /    | |__  |_|\n',
              ' |_____/    |_|    |____|  O')


if __name__ == '__main__':
    m2r = Microbiota2Recon()

    m2r.procedure()
