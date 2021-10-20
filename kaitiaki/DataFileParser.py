class DataFileParser:
    def __init__(self, file_pointer='data'):
        """Parses a datafile. Must be used as a context manager:

        >>> with DataFileParser('data') as dfile:
        >>>    ...

        Keyword Arguments:
            file_pointer {str} -- The location of the datafile (default: {'data'})
        """
        self._file = file_pointer

    def __enter__(self):
        class Parser():
            def __init__(self, file_pointer):
                self._datafile = file_pointer
                self._datafile_pointer = open(self._datafile, 'r+')
                self._original_contents = self._datafile_pointer.read().split("\n")
                self._contents = deepcopy(self._original_contents)

            def _check_scientific_notation(self, param):
                idx = self._get_index_of_parameter(param)
                if idx[0] == 3:
                    return 1 # everything on this line requires scientific notation to 1dp
                elif idx[0] in [17, 18]:
                    return 2 # everything on these lines requires scientific notation to 2dp
                else:
                    if param.lower() in ['ct8', 'ct9', 'ct10']:
                        return 1
                    if param.lower() == 'zs':
                        return 2
                    if param.lower() in ['vrot1', 'vrot2']:
                        return 2
                    if param.lower() in ['facsgmin', 'sgthfac']:
                        return 2
                    if param.lower() in ['hkh', 'gff']:
                        return 2

                return 0

            def _determine_decimal_places(self, param):
                idx = self._get_index_of_parameter(param)
                param = param.lower()
                if idx[0] == 16 and param != 'zs':
                    return 3
                if param in ['fmac', 'fam']:
                    return 2
                if param in ['trc1', 'trc2']:
                    return 1
                if param == 'mwts':
                    return 2
                if idx[0] == 15 and param not in ['ct10', 'ct9', 'ct8']:
                    return 2

                return 0

            def _get_index_of_parameter(self, param):
                """
                                                  .

                                                   .
                                         /^\     .
                                    /\   "V"
                                   /__\   I      O  o
                                  //..\\  I     .S
                                  \].`[/  I
                                  /l\/j\  (]    .  O
                                 /. ~~ ,\/I          .
                                 \\L__j^\/I       o
                                  \/--v}  I     o   .
                                  |    |  I   _________
                                  |    |  I c(`       ')o
                                  |    l  I   \.     ,/
                                _/j  L l\_!  _//^---^\\_
                                     Here be wizard.
                """
                param = param.lower()
                # Structure of the lookup table:
                # 3-tuple: (line number, starting point, ending point)
                lookup_table =  {'nm2':       (0,0,4),      'icl':       (1,0,4),
                                 'it1':       (0,4,8),      'ion':       (1,4,8),
                                 'it2':       (0,8,12),     'iam':       (1,8,12),
                                 'jin':       (0,12,16),    'iop':       (1,12,16),
                                 'out':       (0,16,20),    'inuc':      (1,16,20),
                                 'nch':       (0,20,24),    'ibc':       (1,20,24),
                                 'jp':        (0,24,28),    'icn':       (1,24,28),
                                 'ith':       (0,28,32),    'iml1':      (1,28,32),
                                 'ix':        (0,32,36),    'iml2':      (1,32,36),
                                 'iy':        (0,36,40),    'isgth':     (1,36,40),
                                 'iz':        (0,40,44),    'imo':       (1,40,44),
                                 'imode':     (0,44,None),  'idiff':     (1,44,None),
                                 'nwrt1':     (2,0,4),      'eps':       (3,0,8),
                                 'nwrt2':     (2,4,8),      'del':       (3,8,16),
                                 'nwrt3':     (2,8,12),     'dho':       (3,16,24),
                                 'nwrt4':     (2,12,16),    'dt3':       (3,24,32),
                                 'nwrt5':     (2,16,20),    'ddd':       (3,32,None),

                                 'nsave':     (2,20,24),    'evo:ne1':   (4,0,3),
                                 'nmont':     (2,24,None),  'evo:ne2':   (4,3,6),
                                 'evo:ne3':   (4,6,9),      'evo:nb':    (4,9,12),
                                 'evo:nev':   (4,12,15),    'evo:nf':    (4,15,18),
                                 'evo:j1':    (4,18,21),    'evo:j2':    (4,21,24),
                                 'evo:ih':    (4,24,27),    'evo:jh':    (4,27,None),

                                 'evo:id1':   (5,0,3),      'evo:id2':   (5,3,6),
                                 'evo:id3':   (5,6,9),      'evo:id4':   (5,9,12),
                                 'evo:id5':   (5,12,15),    'evo:id6':   (5,15,18),
                                 'evo:id7':   (5,18,21),    'evo:id8':   (5,21,24),
                                 'evo:id9':   (5,24,27),    'evo:id10':  (5,27,30),
                                 'evo:id11':  (5,30,33),    'evo:id12':  (5,33,36),
                                 'evo:id13':  (5,36,39),    'evo:id14':  (5,39,42),
                                 'evo:id15':  (5,42,45),    'evo:id16':  (5,45,48),
                                 'evo:id17':  (5,48,51),    'evo:id18':  (5,51,54),
                                 'evo:id19':  (5,54,57),    'evo:id20':  (5,57,60),
                                 'evo:id21':  (5,60,63),    'evo:id22':  (5,63,66),
                                 'evo:id23':  (5,66,69),    'evo:id24':  (5,69,72),
                                 'evo:id25':  (5,72,75),    'evo:id26':  (5,75,78),
                                 'evo:id27':  (5,78,81),    'evo:id28':  (5,81,84),
                                 'evo:id29':  (5,84,87),    'evo:id30':  (5,87,None),

                                 'evo:id31':  (6,0,3),      'evo:id32':  (6,3,6),
                                 'evo:id33':  (6,6,9),      'evo:id34':  (6,9,12),
                                 'evo:id35':  (6,12,15),    'evo:id36':  (6,15,18),
                                 'evo:id37':  (6,18,21),    'evo:id38':  (6,21,24),
                                 'evo:id39':  (6,24,27),    'evo:id40':  (6,27,30),
                                 'evo:id41':  (6,30,33),    'evo:id42':  (6,33,36),
                                 'evo:id43':  (6,36,39),    'evo:id44':  (6,39,42),
                                 'evo:id45':  (6,42,45),    'evo:id46':  (6,45,48),
                                 'evo:id47':  (6,48,51),    'evo:id48':  (6,51,54),
                                 'evo:id49':  (6,54,57),    'evo:id50':  (6,57,60),
                                 'evo:id51':  (6,60,63),    'evo:id52':  (6,63,66),
                                 'evo:id53':  (6,66,69),    'evo:id54':  (6,69,72),
                                 'evo:id55':  (6,72,75),    'evo:id56':  (6,75,78),
                                 'evo:id57':  (6,78,81),    'evo:id58':  (6,81,84),
                                 'evo:id59':  (6,84,87),    'evo:id60':  (6,87,None),

                                 'evo:id61':  (7,0,3),      'evo:id62':  (7,3,6),
                                 'evo:id63':  (7,6,9),      'evo:id64':  (7,9,12),
                                 'evo:id65':  (7,12,15),    'evo:id66':  (7,15,18),
                                 'evo:id67':  (7,18,21),    'evo:id68':  (7,21,24),
                                 'evo:id69':  (7,24,27),    'evo:id70':  (7,27,30),
                                 'evo:id71':  (7,30,33),    'evo:id72':  (7,33,36),
                                 'evo:id73':  (7,36,39),    'evo:id74':  (7,39,42),
                                 'evo:id75':  (7,42,45),    'evo:id76':  (7,45,48),
                                 'evo:id77':  (7,48,51),    'evo:id78':  (7,51,54),
                                 'evo:id79':  (7,54,57),    'evo:id80':  (7,57,60),
                                 'evo:id81':  (7,60,63),    'evo:id82':  (7,63,66),
                                 'evo:id83':  (7,66,69),    'evo:id84':  (7,69,72),
                                 'evo:id85':  (7,72,75),    'evo:id86':  (7,75,78),
                                 'evo:id87':  (7,78,81),    'evo:id88':  (7,81,84),
                                 'evo:id89':  (7,84,87),    'evo:id90':  (7,87,None),

                                 'nuc:ne1':   (8,0,3),      'nuc:ne2':   (8,3,6),
                                 'nuc:ne3':   (8,6,9),      'nuc:nb':    (8,9,12),
                                 'nuc:nev':   (8,12,15),    'nuc:nf':    (8,15,18),
                                 'nuc:j1':    (8,18,21),    'nuc:j2':    (8,21,24),
                                 'nuc:ih':    (8,24,27),    'nuc:jh':    (8,27,None),

                                 'nuc:id1':   (9,0,3),      'nuc:id2':   (9,3,6),
                                 'nuc:id3':   (9,6,9),      'nuc:id4':   (9,9,12),
                                 'nuc:id5':   (9,12,15),    'nuc:id6':   (9,15,18),
                                 'nuc:id7':   (9,18,21),    'nuc:id8':   (9,21,24),
                                 'nuc:id9':   (9,24,27),    'nuc:id10':  (9,27,30),
                                 'nuc:id11':  (9,30,33),    'nuc:id12':  (9,33,36),
                                 'nuc:id13':  (9,36,39),    'nuc:id14':  (9,39,42),
                                 'nuc:id15':  (9,42,45),    'nuc:id16':  (9,45,48),
                                 'nuc:id17':  (9,48,51),    'nuc:id18':  (9,51,54),
                                 'nuc:id19':  (9,54,57),    'nuc:id20':  (9,57,60),
                                 'nuc:id21':  (9,60,63),    'nuc:id22':  (9,63,66),
                                 'nuc:id23':  (9,66,69),    'nuc:id24':  (9,69,72),
                                 'nuc:id25':  (9,72,75),    'nuc:id26':  (9,75,78),
                                 'nuc:id27':  (9,78,81),    'nuc:id28':  (9,81,84),
                                 'nuc:id29':  (9,84,87),    'nuc:id30':  (9,87,None),

                                 'nuc:id31':  (10,0,3),     'nuc:id32':  (10,3,6),
                                 'nuc:id33':  (10,6,9),     'nuc:id34':  (10,9,12),
                                 'nuc:id35':  (10,12,15),   'nuc:id36':  (10,15,18),
                                 'nuc:id37':  (10,18,21),   'nuc:id38':  (10,21,24),
                                 'nuc:id39':  (10,24,27),   'nuc:id40':  (10,27,30),
                                 'nuc:id41':  (10,30,33),   'nuc:id42':  (10,33,36),
                                 'nuc:id43':  (10,36,39),   'nuc:id44':  (10,39,42),
                                 'nuc:id45':  (10,42,45),   'nuc:id46':  (10,45,48),
                                 'nuc:id47':  (10,48,51),   'nuc:id48':  (10,51,54),
                                 'nuc:id49':  (10,54,57),   'nuc:id50':  (10,57,60),
                                 'nuc:id51':  (10,60,63),   'nuc:id52':  (10,63,66),
                                 'nuc:id53':  (10,66,69),   'nuc:id54':  (10,69,72),
                                 'nuc:id55':  (10,72,75),   'nuc:id56':  (10,75,78),
                                 'nuc:id57':  (10,78,81),   'nuc:id58':  (10,81,84),
                                 'nuc:id59':  (10,84,87),   'nuc:id60':  (10,87,None),

                                 'nuc:id61':  (11,0,3),     'nuc:id62':  (11,3,6),
                                 'nuc:id63':  (11,6,9),     'nuc:id64':  (11,9,12),
                                 'nuc:id65':  (11,12,15),   'nuc:id66':  (11,15,18),
                                 'nuc:id67':  (11,18,21),   'nuc:id68':  (11,21,24),
                                 'nuc:id69':  (11,24,27),   'nuc:id70':  (11,27,30),
                                 'nuc:id71':  (11,30,33),   'nuc:id72':  (11,33,36),
                                 'nuc:id73':  (11,36,39),   'nuc:id74':  (11,39,42),
                                 'nuc:id75':  (11,42,45),   'nuc:id76':  (11,45,48),
                                 'nuc:id77':  (11,48,51),   'nuc:id78':  (11,51,54),
                                 'nuc:id79':  (11,54,57),   'nuc:id80':  (11,57,60),
                                 'nuc:id81':  (11,60,63),   'nuc:id82':  (11,63,66),
                                 'nuc:id83':  (11,66,69),   'nuc:id84':  (11,69,72),
                                 'nuc:id85':  (11,72,75),   'nuc:id86':  (11,75,78),
                                 'nuc:id87':  (11,78,81),   'nuc:id88':  (11,81,84),
                                 'nuc:id89':  (11,84,87),   'nuc:id90':  (11,87,None),

                                 'isx1':      (12,0,3),     'isx2':      (12,3,6),
                                 'isx3':      (12,6,9),     'isx4':      (12,9,12),
                                 'isx5':      (12,12,15),   'isx6':      (12,15,18),
                                 'isx7':      (12,18,21),   'isx8':      (12,21,24),
                                 'isx9':      (12,24,27),   'isx10':     (12,27,30),
                                 'isx11':     (12,30,33),   'isx12':     (12,33,36),
                                 'isx13':     (12,36,39),   'isx14':     (12,39,42),
                                 'isx15':     (12,42,45),   'isx16':     (12,45,48),
                                 'isx17':     (12,48,51),   'isx18':     (12,51,None),
                                 'isx19':     (13,0,3),     'isx20':     (13,3,6),
                                 'isx21':     (13,6,9),     'isx22':     (13,9,12),
                                 'isx23':     (13,12,15),   'isx24':     (13,15,18),
                                 'isx25':     (13,18,21),   'isx26':     (13,21,24),
                                 'isx27':     (13,24,27),   'isx28':     (13,27,30),
                                 'isx29':     (13,30,33),   'isx30':     (13,33,36),
                                 'isx31':     (13,36,39),   'isx32':     (13,39,42),
                                 'isx33':     (13,42,None),
                                 'isx34':     (14,0,3),     'isx35':     (14,3,6),
                                 'isx36':     (14,6,9),     'isx37':     (14,9,12),
                                 'isx38':     (14,12,15),   'isx39':     (14,15,18),
                                 'isx40':     (14,18,21),   'isx41':     (14,21,24),
                                 'isx42':     (14,24,27),   'isx43':     (14,27,30),
                                 'isx44':     (14,30,33),   'isx45':     (14,33,36),
                                 'isx46':     (14,36,39),   'isx47':     (14,39,42),
                                 'isx48':     (14,42,None),

                                 'dt1':       (15,0,5),     'dt2':       (15,5,10),
                                 'ct1':       (15,10,15),   'ct2':       (15,15,20),
                                 'ct3':       (15,20,25),   'ct4':       (15,25,30),
                                 'ct5':       (15,30,35),   'ct6':       (15,35,40),
                                 'ct7':       (15,40,45),   'ct8':       (15,45,53),
                                 'ct9':       (15,53,61),   'ct10':      (15,61,None),

                                 'zs':        (16,0,9),     'alpha':     (16,9,15),
                                 'ch':        (16,15,21),   'cn':        (16,21,27),
                                 'cn':        (16,27,33),   'co':        (16,33,39),
                                 'cne':       (16,39,45),   'cmg':       (16,45,51),
                                 'csi':       (16,51,57),   'cfe':       (16,57,None),

                                 'rcd':       (17,0,9),     'os':        (17,9,18),
                                 'rml':       (17,18,27),   'rmg':       (17,27,36),
                                 'eca':       (17,36,45),   'xf':        (17,45,54),
                                 'dr':        (17,54,None),

                                 'rmt':       (18,0,9),     'rhl':       (18,9,18),
                                 'ac':        (18,18,27),   'ak1':       (18,27,36),
                                 'ak2':       (18,36,45),   'ect':       (18,45,54),
                                 'trb':       (18,54,None),

                                 'iram':      (19,0,2),     'irs1':      (19,2,4),
                                 'vrot1':     (19,4,13),    'irs2':      (19,13,15),
                                 'vrot2':     (19,15,24),   'fmac':      (19,24,29),
                                 'fam':       (19,29,None),

                                 'ivmc':      (20,0,2),     'trc1':      (20,2,8),
                                 'ivms':      (20,8,10),    'trc2':      (20,10,16),
                                 'mwts':      (20,16,21),   'iagb':      (20,21,23),
                                 'isgfac':    (20,23,25),   'facsgmin':  (20,25,34),
                                 'sgthfac':   (20,34,None),

                                 'istart':    (21,0,2),     'hkh':       (21,2,11),
                                 'gff':       (21,11,None)
                                }

                """
                    God is dead. God remains dead. And we have killed him.

                    - Nietzsche
                """

                if param in lookup_table.keys():
                    return lookup_table[param]
                else:
                    raise KeyError(f"Parameter {param} not recognised.")

            def _write_to_pointer(self, pointer, data, mode):
                if mode == 'replace':
                    pointer.seek(0)
                pointer.write("\n".join(data))

            def make_backup(self):
                with open(self._datafile + ".bak", 'w') as file:
                    self._write_to_pointer(file, self._original_contents)

            def backup_if_not_exists(self):
                if not path.exists(self._datafile + ".bak"):
                    self.make_backup()

            def set(self, param, value):
                idx = self._get_index_of_parameter(param)

                value = str(value)

                num_dp_of_scientific_notation = self._check_scientific_notation(param)
                num_dp = self._determine_decimal_places(param)

                if num_dp_of_scientific_notation > 0:
                    value = f'%.{num_dp_of_scientific_notation}E' % Decimal(value)
                if num_dp > 0:
                    value = '{:.2f}'.format(num_dp)

                if idx[2] == None:
                    endpoint = len(self._contents[idx[0]])
                else:
                    endpoint = idx[2]

                length = endpoint-idx[1]
                value = value.rjust(length)

                leftpart = self._contents[idx[0]][:idx[1]]
                rightpart = self._contents[idx[0]][endpoint:]

                new_line = leftpart + value + rightpart

                self._contents[idx[0]] = new_line

                self._write_to_pointer(self._datafile_pointer, self._contents, 'replace')

            def get(self, param):
                idx = self._get_index_of_parameter(param)

                val = self._contents[idx[0]][idx[1]:idx[2]]

                num_dp = self._check_scientific_notation(param) + self._determine_decimal_places(param)

                if num_dp > 0:
                    return float(val)

                return int(val)

            def restore_backup(self):
                pass

            def set_zams_mass(self, target_mass):
                self.set('RML', target_mass)
                self.set('IML1', 9)

        self._parser = Parser(self._file)

        return self._parser

    def __exit__(self, exc_type, exc_value, traceback):
        self._parser._datafile_pointer.close()