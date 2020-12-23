from pymeeus_oo.calculation.base import TOL, machine_accuracy
from pymeeus_oo.calculation.base import get_ordinal_suffix


def main():

    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's print the tolerance
    print_me("The default value for the tolerance is", TOL)

    # Find the accuracy of this computer
    j, d = machine_accuracy()
    print_me("Number of significant BITS in the mantissa\t", j)
    print_me("Number of significant DIGITS in a decimal number", d)

    print("")

    print_me("The suffix for ordinal 2 is", get_ordinal_suffix(2))
    print_me("The suffix for ordinal 11 is", get_ordinal_suffix(11))
    print_me("The suffix for ordinal 12 is", get_ordinal_suffix(12))
    print_me("The suffix for ordinal 13 is", get_ordinal_suffix(13))
    print_me("The suffix for ordinal 14 is", get_ordinal_suffix(14))
    print_me("The suffix for ordinal 16 is", get_ordinal_suffix(16))
    print_me("The suffix for ordinal 23 is", get_ordinal_suffix(23))


if __name__ == "__main__":
    main()
