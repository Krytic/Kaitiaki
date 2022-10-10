import kaitiaki

COMMENT_CHARACTER = '#'
BEGIN_BLOCK_CHAR = '['
END_BLOCK_CHAR = ']'
KEY_VALUE_SEPARATOR = '='
MULTILINE_COMMENT = '"""' # The end character is MULTILINE_COMMENT[::-1]
START_TOKEN = "begin"
END_TOKEN = "end"

class Lexer:
    def __init__(self, filepath):
        self.load_file_to_dictionary(filepath)
        validation = self.validate_lexicon()

        if not validation[0]:
            raise ValueError(f"Invalid structure to lexicon (\"{validation[1][1]}\" is not a {validation[1][0]} key).")

    def __str__(self):
        plural = ""
        num = len(self._data.keys())
        if num > 1: plural = "s"
        return f"<<kaitiaki lexicon object representing {num} run{plural}>>"

    def load_file_to_dictionary(self, filepath):
        with open(filepath, 'r') as f:
            block = -1

            data = dict()
            data['options'] = kaitiaki.constants.DEFAULT_LEXICON_OPTIONS

            in_multiline_comment = False
            in_option_block = False

            for line in f.readlines():
                line = line.strip()

                if len(line) == 0: continue
                if line[0:len(MULTILINE_COMMENT)] == MULTILINE_COMMENT and not in_multiline_comment:
                    in_multiline_comment = True
                    continue

                if in_multiline_comment:
                    if line[0:len(MULTILINE_COMMENT)] == MULTILINE_COMMENT[::-1]:
                        in_multiline_comment = False
                    continue

                if line[0] == COMMENT_CHARACTER: continue

                if line.strip() == f"{BEGIN_BLOCK_CHAR}options{END_BLOCK_CHAR}":
                    in_option_block = True
                elif line.strip() == f"{BEGIN_BLOCK_CHAR}end options{END_BLOCK_CHAR}":
                    in_option_block = False
                elif line[0] == BEGIN_BLOCK_CHAR:
                    endblock = line.find(END_BLOCK_CHAR)

                    blockname = line[1:endblock].strip()
                    block += 1

                    current_block_name = blockname

                    data[block] = dict()
                else:
                    parts = line.split(KEY_VALUE_SEPARATOR)

                    parameter, value = parts[0].strip(), parts[1].strip()

                    key = block
                    if in_option_block: key = 'options'
                    if in_option_block and parameter == 'is_binary':
                        parameter = bool(parameter)

                    data[key][parameter] = value

        self._data = data

    def validate_lexicon(self):
        for block, data in self._data.items():
            if block == 'options':
                # do option block validation
                continue
            else:
                for key, value in data.items():
                    if key.lower() not in kaitiaki.constants.dfile_struct.keys():
                        return False, ('dfile', key)
        return True, None

    def fetch_lexicon(self):
        return self._data

    def __iter__(self):
        self.__cntr = 0
        return self

    def __next__(self):
        if self.__cntr >= len(self._data):
            raise StopIteration
        else:
            data = self._data[self.__cntr]
            self.__cntr += 1

            return data