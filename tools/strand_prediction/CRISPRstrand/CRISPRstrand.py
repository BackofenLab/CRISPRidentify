#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import preprocessing as pp
import utils as u

if __name__ == "__main__":
    cmdline_parser = argparse.ArgumentParser('CRISPRstrand_v2')

    cmdline_parser.add_argument('-cv', '--cross_validation',
                                action='store_true',
                                default=False,
                                help='Cross validation (only applied for training)')

    cmdline_parser.add_argument('-tr', '--training',
                                action='store_true',
                                default=False,
                                help='Whether to train a model')

    cmdline_parser.add_argument('-i', '--input_files',
                                nargs='+', 
                                default= ['Example/Input3.fa'],
                                help='Filenames of the input data.')

    cmdline_parser.add_argument('-cols', '--usecols',
                                nargs='+', 
                                default= [0, 5],
                                help='ID and consensus repeat fields to use. Must be specified for csv/tsv/xls...')

    cmdline_parser.add_argument('-m', '--model_path',
                                default='Models/model_r.h5',
                                help='Evaluation/prediction model path',
                                type=str)

    cmdline_parser.add_argument('-type', '--repeat_type',
                                action='store_true',
                                default=False,
                                help='Whether to train a model')

    cmdline_parser.add_argument('-out', '--output_folder',
                                default='Results',
                                help='Output save',
                                type=str)


    args, unknowns = cmdline_parser.parse_known_args()
    print(args)

    tr = args.training
    cv = args.cross_validation
    inputs = args.input_files
    columns = None if len(args.usecols) == 0 else args.usecols
    model_path = args.model_path
    do_type = args.repeat_type
    output_folder = args.output_folder
    ###################################
    
    
    if tr:
        
        from keras.backend import clear_session
        import train.train as tr

        df = pp.prepare_input(inputs, usecols = columns)
        X, y, sequence_length = pp.process_dataset_for_training(df)
            
        run_dict = {'batch_size'        : 128,
                    'num_epochs'        : 100,
                    'num_repeats'       : 1,
                    'k'                 : 5, # only will be used in CV
                    'classifier_build'  : 'parallel',
                    'kernel_width'      : [4, 5, 6, 8, 10], # put 6
                    'sequence_length'   : sequence_length
                    }

        clear_session()

        if cv:
            classifiers, histories, roc_auc_scores = te.cv_run(X, y, run_dict)
        else:
            classifiers, histories, roc_auc_scores = te.run(X, y, run_dict)

    else:

        import evaluate as ev

        classifier, sequence_length = u.load_model(model_path)
        df = pp.prepare_input(inputs, output_folder, False, usecols = columns, do_type = do_type)
        df, X, y = pp.process_dataset_for_test(df, sequence_length)
        df_out, roc_auc_scores = ev.test(classifier, df, X, y, output_folder)

    print(roc_auc_scores)
