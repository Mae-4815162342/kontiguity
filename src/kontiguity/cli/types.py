import click

import cooler
import os

# # accepts a list of cool files or a file containing a list of cool, one per line
# class coolType(click.ParamType):
#     name="cool"

#     def is_cool(self, file):
#         try:
#             return cooler.fileops.is_cooler(file)
#         except:
#             return False
        
#     def is_mcool(self, file):
#         try:
#             return cooler.fileops.is_multires_file(file)
#         except:
#             return False

#     def convert(self, value, param, ctx):
#         if self.is_cool(value) or self.is_mcool(value):
#             return [value]
#         if not os.path.isfile(value):
#             files = value.split(',')
#             for file in files:
#                 if not self.is_cool(file) and not self.is_mcool(file):
#                     self.fail(f"{file} is not a valid cool file", param, ctx)
#             return files
#         else:
#             files = []
#             with open(value, 'r') as f:
#                 for line in f.readlines():
#                     file = line.replace('\n', '').replace(' ', '')
#                     if not self.is_cool(file) and not self.is_mcool(file):
#                         self.fail(f"{file} is not a valid cool file", param, ctx)
#                     else:
#                         files.append(file)
#             return files

# class IntListType(click.ParamType):
#     name="int_list"

#     def convert(self, value, param, ctx):
#         list = value.split(',')
#         ints = []
#         for i in list:
#             if not isinstance(int(i), int):
#                 self.fail(f"{value} is not a comma-separated integer list.")
#             if int(i) < 0:
#                 self.fail(f"All integers provided in {param.name} must be positive integers.")
#             ints.append(int(i))
#         return ints
    
# class StrListType(click.ParamType):
#     name="str_list"

#     def convert(self, value, param, ctx):
#         strs = value.split(',')
#         for s in strs:
#             if len(s) == 0:
#                 self.fail(f"Empty strings. Comma must separate two distinct values; if a single string is passed, no comma must be found in {param.name}.")
#         return strs

    
class PairListType(click.ParamType):
    name="pair_list"

    def convert(self, value, param, ctx):
        strs = value.split(',')
        pairs = []
        for s in strs:
            if len(s) == 0:
                self.fail(f"Empty strings. Comma must separate at least two distinct values; if a single string is passed, no comma must be found in {param.name}.")
            else:
                pair = s.split(':')
                if len(pair) > 2 or len(pair) == 0:
                    self.fail(f"Wrong number of values provided in {pair}. A single value or a colon separated pair was expected.")
                elif len(pair) == 2:
                    if len(pair[0]) == 0 or len(pair[1]) == 0:
                        self.fail(f"Empty strings. Colon must separate at least two distinct values; if a single value is passed, no colon must be found in {param.name}.")
                    else:
                        pairs.append((pair[0], pair[1], "paired"))
                else:
                    if len(pair[0]) == 0:
                        self.fail(f"Empty strings. Colon must separate at least two distinct values; if a single value is passed, no colon must be found in {param.name}.")
                    else:
                        pairs.append((pair[0],"", "unpaired"))
            
        return pairs

# COOL = coolType()
# INT_LIST = IntListType()
# STR_LIST = StrListType()
PAIR_LIST = PairListType()