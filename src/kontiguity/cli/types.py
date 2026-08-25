import click

class IntListType(click.ParamType):
    name="int_list"

    def convert(self, value, param, ctx):
        list = value.split(',')
        ints = []
        for i in list:
            if not isinstance(int(i), int):
                self.fail(f"{value} is not a comma-separated integer list.")
            if int(i) < 0:
                self.fail(f"All integers provided in {param.name} must be positive integers.")
            ints.append(int(i))
        return ints
    def convert(self, value, param, ctx):
        list = value.split(',')
        ints = []
        for i in list:
            if not isinstance(int(i), int):
                self.fail(f"{value} is not a comma-separated integer list.")
            if int(i) < 0:
                self.fail(f"All integers provided in {param.name} must be positive integers.")
            ints.append(int(i))
        return ints
    
class StrListType(click.ParamType):
    name="str_list"

    def convert(self, value, param, ctx):
        strs = value.split(',')
        for s in strs:
            if len(s) == 0:
                self.fail(f"Empty strings. Comma must separate two distinct values; if a single string is passed, no comma must be found in {param.name}.")
        return strs
    
class PairListType(click.ParamType):
    name="pair_list"

    def convert(self, value, param, ctx):
        pairs = value.split(',')
        formated_pairs = []
        for pair in pairs:
            if len(pair) == 0:
                self.fail(f"Empty strings. Comma must separate two distinct values; if a single string is passed, no comma must be found in {param.name}.")
            pair_values = pair.split(':')
            if len(pair_values) > 2:
                self.fail(f"Too many values provided in {pair}. Colons must separate at most two distinct values; if a single value is passed, do not add colons.")
            is_paired = len(pair_values) == 2
            formated_pairs.append((pair_values, is_paired))
        return formated_pairs

INT_LIST = IntListType()
STR_LIST = StrListType()
PAIR_LIST = PairListType()