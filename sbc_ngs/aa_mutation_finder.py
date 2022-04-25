
import sys
import csv
import os
from Bio import SeqIO
from Bio.Seq import MutableSeq
import pandas as pd
import ast
import re
from datetime import datetime
from typing import Tuple, List

class MutationFinder():

    def __init__(self, output_path: str, result_folder_id: str, ref_seq_fasta_name: str):

        self.__output_results_folder_path = os.path.join(output_path, "results", result_folder_id)
        self.__ref_seq = MutableSeq(self.read_ref_seq(output_path, ref_seq_fasta_name))
        self.__summary_df = self.read_raw_summary()

        # Add columns to summary_df
        empty_list = ["[]"] * len(self.__summary_df)

        self.__summary_df.insert(12, "aa_mutations", empty_list, True)
        self.__summary_df.insert(13, "aa_nucleotides", empty_list, True)
        self.__summary_df.insert(14, "aa_indels", empty_list, True)
        self.__summary_df.insert(15, "aa_deletions", empty_list, True)

    def read_mutation(self, mutation: str) -> Tuple[str, int, str]:
        ''' Read parameters from mutation in summary e.g. A23G, GTA45AGG'''

        match = re.search('([CGATN*]+)(\d+)([CGATN*]+)', str(mutation))

        return match.group(1), int(match.group(2)), match.group(3)

    def process_summary(self):

        #Iterate through each row.
        for row in range(len(self.__summary_df)):
            print(self.__summary_df['forward'][row] + "_" + self.__summary_df['reverse'][row] + "_" + self.__summary_df['barcode_type'][row])
        # Find mutants from 
        # a. mutation column
            print("Mutations " + datetime.now().strftime("%H:%M:%S"))
            results = self.convert_entry_mutations('mutations', row)
            self.__summary_df['aa_mutations'][row] = results

        # b. nucleotides column 
            print("Nucleotides " + datetime.now().strftime("%H:%M:%S"))
            results = self.convert_entry_mutations('nucleotides', row)
            self.__summary_df['aa_nucleotides'][row] = results

        # c. indels column
            print("InDels " + datetime.now().strftime("%H:%M:%S"))
            results = self.convert_entry_mutations('indels', row)
            self.__summary_df['aa_indels'][row] = results

        # d. deletions column
            print("Deletions " + datetime.now().strftime("%H:%M:%S"))
            results = self.convert_entry_deletions(row)
            self.__summary_df['aa_deletions'][row] = results

        self.generate_output()

    def convert_entry_deletions(self, row) -> List[List[str]]:
        results = []
        if self.__summary_df['deletions'][row]:
            deletions = list(ast.literal_eval(self.__summary_df['deletions'][row]))
            for deletion in deletions:
                # Convert the string or single value to list to list 
                match = re.search('(\d+)-(\d+)', str(deletion))
                if match:
                    deletion_start = int(match.group(1))
                    deletion_end = int(match.group(2))
                else:
                    deletion_start = int(deletion)
                    deletion_end = int(deletion)

                mut_seq = self.generate_deletion_sequence(deletion_start, deletion_end)

                mutation_result = self.find_aa_differences(mut_seq)

                # To handle the potential for multiple results
                if mutation_result:
                    results.append(mutation_result)

        return results

    def convert_entry_mutations(self, column: str, row: int) -> List[Tuple[List[str],str]]:
        results = []
        if self.__summary_df[column][row]:
            mutations = list(ast.literal_eval(self.__summary_df[column][row]))
            
            for mutation in mutations:
                mut_seq = self.generate_mutant(mutation[0])
                mutation_result = self.find_aa_differences(mut_seq)

                # To handle the potential for multiple results
                if mutation_result:
                    results.append((mutation_result, mutation[1]))

        return results

    def generate_deletion_sequence(self, deletion_start: int, deletion_end: int) -> str:
        '''Generate deletion sequence from mismatch.'''

        # Create mutant protein sequence
        mut_seq = MutableSeq(self.__ref_seq)

        pre_seq = mut_seq[:deletion_start-1]
        post_seq = mut_seq[deletion_end:]

        return pre_seq + post_seq

    def generate_mutant(self, mutation: str) -> str: 
        '''Generate mutant sequence from mismatch.'''
        cur, pos, alt = self.read_mutation(mutation)

        cur_len = len(cur)
        alt_len = len(alt)

        # Create mutant protein sequence
        mut_seq = MutableSeq(self.__ref_seq)

        pre_seq = mut_seq[:pos-1]
        site = mut_seq[pos-1:pos-1+cur_len]
        post_seq = mut_seq[pos+cur_len:]

        #Need to confirm it matches original

        if site != cur:
            print("Unmatching site, expecting " + cur + " got " + site)
        else:
            mut_seq = pre_seq + alt + post_seq

        return mut_seq

    def find_aa_differences(self, mut_seq: str) -> List[str]:
        '''Translate mutant and reference sequence and compare them directly. Returns list of identified mismatches.'''

        mismatch_list = []

        # NEED TO MANAGE PARTIAL CODONS

        #if (len(ref_seq) % 3) == 1:
        #    ref_seq.append('N')
        #    ref_seq.append('N')
        #elif (len(ref_seq) % 3) == 2:
        #    ref_seq.append('N')

        ref_aa_seq = self.__ref_seq.translate()
        mut_aa_seq = mut_seq.translate()
        #print(ref_seq.translate(cds=True))

        # Trim to match a.a. lengths
        if len(ref_aa_seq) > len(mut_aa_seq):
            ref_aa_seq = ref_aa_seq[:len(mut_aa_seq)]
        elif len(ref_aa_seq) < len(mut_aa_seq):
            mut_aa_seq = mut_aa_seq[:len(ref_aa_seq)]

        cur_ref_aa_seq = ""
        cur_mut_aa_seq = ""
        starting_position = 0

        for current_position in range(len(ref_aa_seq)):
            if ref_aa_seq[current_position] != mut_aa_seq[current_position]:
                cur_ref_aa_seq = cur_ref_aa_seq + ref_aa_seq[current_position]
                cur_mut_aa_seq = cur_mut_aa_seq + mut_aa_seq[current_position]
            else:
                # if sequence not empty save off    
                if cur_ref_aa_seq != "" or cur_mut_aa_seq != "":
                    mutant = cur_ref_aa_seq + str(starting_position) + cur_mut_aa_seq
                    mismatch_list.append(mutant)
                    cur_ref_aa_seq = ""
                    cur_mut_aa_seq = ""

                # Reset starting position to next.
                starting_position = current_position + 1

        # Need to handle if mutation stored when loop ended.
        if cur_ref_aa_seq != "" or cur_mut_aa_seq != "":
            mutant = cur_ref_aa_seq + str(starting_position) + cur_mut_aa_seq
            mismatch_list.append(mutant)

        return mismatch_list

    def generate_output(self):

        #for row in range(0, 20):
        #    for col in ['aa_nucleotides','aa_indels','aa_deletions']:
        #        entry = list(ast.literal_eval(self.__summary_df[col][row]))
                #length = len(entry)
        #        print(str(row) + col + " " + str(len(entry)) + " " + entry)

        '''Output updated dataframe as csv.'''
        raw_aa_summary_path = os.path.join(self.__output_results_folder_path, "raw_aa_summary.csv")
        self.__summary_df.to_csv(raw_aa_summary_path)

    def read_raw_summary(self) -> pd.DataFrame:
        raw_summary_path = os.path.join(self.__output_results_folder_path, "raw_summary.csv")
        df = pd.read_csv(raw_summary_path)

        return df
    
    def read_ref_seq(self, output_path: str, ref_seq_file_name: str) -> str:

        ref_seq_path = os.path.join(output_path, 'seqs', ref_seq_file_name)

        for record in SeqIO.parse(ref_seq_path, "fasta"):
            seq = record.seq

        return seq

def main(args):
    '''main method'''
    finder = MutationFinder(output_path=args[0], result_folder_id=args[1], ref_seq_fasta_name=args[2])
    finder.process_summary()

if __name__ == '__main__':

    output_path = "/home/will/Source/SequenceGenie/example/fasta/"
    result_folder_id = "9d81e687-6952-4042-9c7f-5d6812e3205c"
    ref_seq_fasta_name = "SBC003382.fasta"

    main([output_path, result_folder_id, ref_seq_fasta_name])
    # main(sys.argv[1:])





