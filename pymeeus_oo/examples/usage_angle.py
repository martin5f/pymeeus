from math import (
    pi
)

from pymeeus_oo.calculation.angle import Angle


def main():

    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's show some uses of Angle class
    print("\n" + 35 * "*")
    print("*** Use of Angle class")
    print(35 * "*" + "\n")

    # Create an Angle object, providing degrees, minutes and seconds
    a = Angle(-23.0, 26.0, 48.999983999)
    print("a = Angle(-23.0, 26.0, 48.999983999)")

    # First we print using the __call__ method (note the extra parentheses)
    print_me("The angle 'a()' is", a())  # -23.44694444

    # Second we print using the __str__ method (no extra parentheses needed)
    print_me("The angle 'a' is", a)  # -23.44694444

    print("")

    # Use the copy constructor
    b = Angle(a)
    print_me("Angle 'b', which is a copy of 'a', is", b)

    print("")

    # Use the static 'deg2dms()' method to carry out conversions
    d, m, s, sign = Angle.deg2dms(23.44694444)
    val = "{}d {}' {}''".format(int(sign * d), m, s)
    print_me("{Deg}d {Min}' {Sec}''", val)  # 23d 26' 48.999984''

    # We can print Angle 'a' directly in sexagesimal format
    # In 'fancy' format:                                # -23d 26' 48.999984''
    print_me("{Deg}d {Min}' {Sec}''", a.dms_str(n_dec=6))
    # In plain format:
    print_me("{Deg}:{Min}:{Sec}", a.dms_str(False, 6))  # -23:26:48.999983999

    print("")

    # Print directly as a tuple
    a = Angle(23.44694444)
    print_me("a.dms_tuple()", a.dms_tuple())
    print_me("a.ra_tuple()", a.ra_tuple())

    print("")

    # Redefine Angle 'a' several times
    a.set(-0.44694444)
    print("a.set(-0.44694444)")
    print_me("   a.dms_str()", a.dms_str())  # -26' 48.999984''
    a.set(0, 0, -46.31)
    print("a.set(0, 0, -46.31)")
    print_me("   a.dms_str(False)", a.dms_str(False))  # 0:0:-46.31

    print("")

    # We can use decimals in degrees/minutes. They are converted automatically
    a.set(0, -46.25, 0.0)
    print("a.set(0, -46.25, 0.0)")
    print_me("   a.dms_str()", a.dms_str())  # -46' 15.0''
    a.set(0, 0, 0.0)
    print("a.set(0, 0, 0.0)")
    print_me("   a.dms_str()", a.dms_str())  # 0d 0' 0.0''

    print("")

    # We can define the angle as in radians. It will be converted to degrees
    b = Angle(pi, radians=True)
    print_me("b = Angle(pi, radians=True); print(b)", b)  # 180.0

    # And we can easily carry out the 'degrees to radians' conversion
    print_me("print(b.rad())", b.rad())  # 3.14159265359

    print("")

    # We can also specify the angle as a Right Ascension
    print("Angle can be given as a Right Ascension: Hours, Minutes, Seconds")
    a.set_ra(9, 14, 55.8)
    print("a.set_ra(9, 14, 55.8)")
    print_me("   print(a)", a)
    b = Angle(9, 14, 55.8, ra=True)
    print("b = Angle(9, 14, 55.8, ra=True)")
    print_me("   print(b)", b)

    print("")

    # We can print the Angle as Right Ascension, as a float and as string
    a = Angle(138.75)
    print("a = Angle(138.75)")
    print_me("   print(a.get_ra())", a.get_ra())
    print_me("   print(a.ra_str())", a.ra_str())
    print_me("   print(a.ra_str(False))", a.ra_str(False))

    print("")

    # Use the 'to_positive()' method to get the positive version of an angle
    a = Angle(-87.32)  # 272.68
    print("a = Angle(-87.32)")
    print_me("   print(a.to_positive())", a.to_positive())

    print("")

    # Call the __repr__() method to get a string defining the current object
    # This string can then be fed to 'eval()' function to generate the object
    print_me("print(b.__repr__())", b.__repr__())  # Angle(138.7325)
    c = eval(repr(b))
    print_me("c = eval(repr(b)); print(c)", c)  # 138.7325

    print("")

    print_me("c", c)  # 138.7325

    # Negate an angle
    d = Angle(13, 30)
    print_me("d", d)  # 13.5
    e = -d
    print_me("   e = -d", e)  # -13.5

    # Get the absolute value of an angle
    e = abs(e)
    print_me("   e = abs(e)", e)  # 13.5

    # Module operation on an angle
    d = Angle(17.0)
    print_me("d", d)  # 17.0
    e = c % d
    print_me("   e = c % d", e)  # 10.0

    print("")

    # Convert the angle to an integer
    d = Angle(13.95)
    print_me("d", d)  # 13.95
    print_me("   int(d)", int(d))  # 13.0
    d = Angle(-4.95)
    print_me("d", d)  # -4.95
    print_me("   int(d)", int(d))  # -4.0

    # Convert the angle to a float
    print_me("   float(d)", float(d))  # -4.95

    # Round the angle to a float
    e = Angle(-4.951648)
    print_me("e", e)  # -4.951648
    print_me("   round(e)", round(e))  # -5.0
    print_me("   round(e, 2)", round(e, 2))  # -4.95
    print_me("   round(e, 3)", round(e, 3))  # -4.952
    print_me("   round(e, 4)", round(e, 4))  # -4.9516

    print("")

    # Comparison operators
    print_me("   d == e", d == e)  # False
    print_me("   d != e", d != e)  # True
    print_me("   d > e", d > e)  # True
    print_me("   c >= 180.0", c >= 180.0)  # False
    print_me("   c < 180.0", c < 180.0)  # True
    print_me("   c <= 180.0", c <= 180.0)  # True

    print("")

    # It is very easy to add Angles to obtain a new Angle
    e = c + d
    print_me("   c + d", e)  # 133.7825

    # We can also directly add a decimal angle
    e = c + 11.5
    print_me("   c + 11.5", e)  # 150.2325

    print("")

    # Types allowed are int, float and Angle
    print('e = c + "32.5"')
    try:
        e = c + "32.5"
    except TypeError:
        print("TypeError!: Valid types are int, float, and Angle, not string!")

    print("")

    # Subtraction
    e = c - d
    print_me("   c - d", e)  # 143.6825

    # Multiplication
    c.set(150.0)
    d.set(5.0)
    print_me("c", c)  # 150.0
    print_me("d", d)  # 5.0
    e = c * d
    print_me("   c * d", e)  # 30.0

    # Division
    c.set(150.0)
    d.set(6.0)
    print_me("d", d)  # 6.0
    e = c / d
    print_me("   c / d", e)  # 25.0

    print("")

    # Division by zero is not allowed
    d.set(0.0)
    print_me("d", d)  # 0.0
    print("e = c / d")
    try:
        e = c / d
    except ZeroDivisionError:
        print("ZeroDivisionError!: Division by zero is not allowed!")

    print("")

    # Power
    d.set(2.2)
    print_me("d", d)  # 2.2
    e = c ** d
    print_me("   c ** d", e)  # 91.57336709992524

    print("")

    # Accumulative module operation
    d.set(17.0)
    print_me("d", d)  # 17.0
    e %= d
    print_me("   e %= d", e)  # 6.573367099925235

    # Accumulative addition
    c += d
    print_me("   c += d", c)  # 167.0

    # Accumulative subtraction
    print_me("b", b)  # 138.7325
    c -= b
    print_me("   c -= b", c)  # 28.2675

    # Accumulative multiplication
    print_me("b", b)  # 138.7325
    c *= b
    print_me("   c *= b", c)  # 321.62094375

    # Accumulative division
    print_me("b", b)  # 138.7325
    d.set(6.0)
    print_me("d", d)  # 6.0
    b /= d
    print_me("   b /= d", b)  # 23.1220833333

    # Accumulative power
    d.set(2.2)
    print_me("d", d)  # 2.2
    c = abs(c)
    print_me("   c = abs(c)", c)  # 321.62094375
    c **= d
    print_me("   c **= d", c)  # 254.307104203

    print("")

    # The same operation, but by the right side
    e = 3.5 + b
    print_me("   e = 3.5 + b", e)  # 26.6220833333
    e = 3.5 - b
    print_me("   e = 3.5 - b", e)  # -19.6220833333
    e = 3.5 * b
    print_me("   e = 3.5 * b", e)  # 80.9272916667
    e = 3.5 / b
    print_me("   e = 3.5 / b", e)  # 0.151370443119
    e = 3.5 ** b
    print_me("   e = 3.5 ** b", e)  # 260.783691406


if __name__ == "__main__":
    main()
