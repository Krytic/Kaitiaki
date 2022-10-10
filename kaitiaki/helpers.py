import kaitiaki

class Range:
    def __init__(self, bound_low, bound_high, inclusivity="[]", mutable=False):
        self._bound_low = bound_low
        self._bound_high = bound_high

        self._right_inclusive = inclusivity[1] == "]"

        self._inclusives = [None, None]

        if inclusivity[0] == "[":
            self._left_inclusive = True
            self._inclusives[0] = "["
        else:
            self._left_inclusive = False
            self._inclusives[0] = "("

        if inclusivity[1] == "]":
            self._right_inclusive = True
            self._inclusives[1] = "]"
        else:
            self._right_inclusive = False
            self._inclusives[1] = ")"

        self._mutable = mutable

    def __contains__(self, item):
        part1 = False
        part2 = False

        if self._left_inclusive:
            part1 = self._bound_low <= item
        else:
            part1 = self._bound_low <  item

        if self._right_inclusive:
            part2 = self._bound_high >= item
        else:
            part2 = self._bound_high >  item

        return (part1 and part2)

    def __getitem__(self, idx):
        if idx == 0:
            return self._bound_low
        elif idx == 1:
            return self._bound_high

        raise ValueError("idx not allowed -- use 0 to access lower bound, or 1 to access upper bound")

    def inclusives(self):
        return (self._left_inclusive, self._right_inclusive)

    def __setitem__(self, key, val):
        if self._mutable:
            if key == 0:
                self._bound_low = val
            elif key == 1:
                self._bound_high = val
            else:
                raise ValueError("Key was set improperly. Use 0 to modify the lower bound, or 1 to modify the upper bound.")
        else:
            raise TypeError("This object was not initialized with mutable methods.")

    def __str__(self):
        return f"{self._inclusives[0]}{self[0]}, {self[1]}{self._inclusives[1]}"