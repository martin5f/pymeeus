# -*- coding: utf-8 -*-


# PyPlanets: Object-oriented refactoring of PyMeeus, a Python library implementing astronomical algorithms.
# Copyright (C) 2020  Martin Fünffinger
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


import datetime

from pyplanets.core.angle import Angle
from pyplanets.core.base import get_ordinal_suffix, iint
from pyplanets.core.epoch import DAY2MIN
from pyplanets.core.epoch import Epoch


def main():

    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's do some work with the Epoch class
    print("\n" + 35 * "*")
    print("*** Use of Epoch class")
    print(35 * "*" + "\n")

    e = Epoch(1987, 6, 19.5)
    print_me("JDE for 1987/6/19.5", e)

    # Redefine the Epoch object
    e.set(333, "Jan", 27, 12)
    print_me("JDE for 333/1/27.5", e)

    # We can create an Epoch from a 'date' or 'datetime' object
    d = datetime.datetime(837, 4, 10, 7, 12, 0, 0)
    f = Epoch(d)
    print_me("JDE for 837/4/10.3", f)

    print("")

    # Check if a given date belong to the Julian or Gregorian calendar
    print_me("Is 1590/4/21.4 a Julian date?", Epoch.is_julian(1590, 4, 21.4))

    print("")

    # We can also check if a given year is leap or not
    print_me("Is -1000 a leap year?", Epoch.is_leap(-1000))
    print_me("Is 1800 a leap year?", Epoch.is_leap(1800))
    print_me("Is 2012 a leap year?", Epoch.is_leap(2012))

    print("")

    # Get the Day Of Year corresponding to a given date
    print_me("Day Of Year (DOY) of 1978/11/14", Epoch.get_doy(1978, 11, 14))
    print_me("Day Of Year (DOY) of -400/2/29.9", Epoch.get_doy(-400, 2, 29.9))

    print("")

    # Now the opposite: Get a date from a DOY
    t = Epoch.doy2date(2017, 365.7)
    s = str(t[0]) + "/" + str(t[1]) + "/" + str(round(t[2], 2))
    print_me("Date from DOY 2017:365.7", s)

    t = Epoch.doy2date(-4, 60)
    s = str(t[0]) + "/" + str(t[1]) + "/" + str(round(t[2], 2))
    print_me("Date from DOY -4:60", s)

    print("")

    # There is an internal table which we can use to get the leap seconds
    print_me("Number of leap seconds applied up to July 1983",
             Epoch.leap_seconds(1983, 7))

    print("")

    # We can convert the internal JDE value back to a date
    e = Epoch(2436116.31)
    y, m, d = e.get_date()
    s = str(y) + "/" + str(m) + "/" + str(round(d, 2))
    print_me("Date from JDE 2436116.31", s)

    print("")

    # It is possible to get the day of the week corresponding to a given date
    e = Epoch(2018, "Feb", 15)
    print_me("The day of week of 2018/2/15 is", e.dow(as_string=True))

    print("")

    # In some cases it is useful to get the Modified Julian Day (MJD)
    e = Epoch(1923, "August", 23)
    print_me("Modified Julian Day for 1923/8/23", round(e.mjd(), 2))

    print("")

    # If your system is appropriately configured, you can get the difference in
    # seconds between your local time and UTC
    print_me(
        "To convert from local system time to UTC you must add/subtract"
        + " this amount of seconds",
        Epoch.utc2local(),
    )

    print("")

    # Compute DeltaT = TT - UT differences for various dates
    print_me("DeltaT (TT - UT) for Feb/333", round(Epoch.tt2ut(333, 2), 1))
    print_me("DeltaT (TT - UT) for Jan/1642", round(Epoch.tt2ut(1642, 1), 1))
    print_me("DeltaT (TT - UT) for Feb/1928", round(Epoch.tt2ut(1928, 1), 1))
    print_me("DeltaT (TT - UT) for Feb/1977", round(Epoch.tt2ut(1977, 2), 1))
    print_me("DeltaT (TT - UT) for Jan/1998", round(Epoch.tt2ut(1998, 1), 1))

    print("")

    # The difference between civil day and sidereal day is almost 4 minutes
    e = Epoch(1987, 4, 10)
    st1 = round(e.mean_sidereal_time(), 9)
    e = Epoch(1987, 4, 11)
    st2 = round(e.mean_sidereal_time(), 9)
    ds = (st2 - st1) * DAY2MIN
    msg = "{}m {}s".format(iint(ds), (ds % 1) * 60.0)
    print_me("Difference between sidereal time 1987/4/11 and 1987/4/10", msg)

    print("")

    print(
        "When correcting for nutation-related effects, we get the "
        + "'apparent' sidereal time:"
    )
    e = Epoch(1987, 4, 10)
    print("e = Epoch(1987, 4, 10)")
    print_me(
        "e.apparent_sidereal_time(23.44357, (-3.788)/3600.0)",
        e.apparent_sidereal_time(23.44357, (-3.788) / 3600.0),
    )
    #    0.549145082637

    print("")

    # Epoch class can also provide the date of Easter for a given year
    # Let's spice up the output a little bit, calling dow() and get_month()
    month, day = Epoch.easter(2019)
    e = Epoch(2019, month, day)
    s = (
        e.dow(as_string=True)
        + ", "
        + str(day)
        + get_ordinal_suffix(day)
        + " of "
        + Epoch.get_month(month, as_string=True)
    )
    print_me("Easter day for 2019", s)
    # I know Easter is always on Sunday, by the way... ;-)

    print("")

    # Compute the date of the Jewish Easter (Pesach) for a given year
    month, day = Epoch.jewish_pesach(1990)
    s = (
        str(day)
        + get_ordinal_suffix(day)
        + " of "
        + Epoch.get_month(month, as_string=True)
    )
    print_me("Jewish Pesach day for 1990", s)

    print("")

    # Now, convert a date in the Moslem calendar to the Gregorian calendar
    y, m, d = Epoch.moslem2gregorian(1421, 1, 1)
    print_me("The date 1421/1/1 in the Moslem calendar is, in Gregorian " +
             "calendar", "{}/{}/{}".format(y, m, d))
    y, m, d = Epoch.moslem2gregorian(1439, 9, 1)
    print_me(
        "The start of Ramadan month (9/1) for Gregorian year 2018 is",
        "{}/{}/{}".format(y, m, d),
    )
    # We can go from the Gregorian calendar back to the Moslem calendar too
    print_me(
        "Date 1991/8/13 in Gregorian calendar is, in Moslem calendar",
        "{}/{}/{}".format(*Epoch.gregorian2moslem(1991, 8, 13)),
    )
    # Note: The '*' before 'Epoch' will _unpack_ the tuple into components

    print("")

    # It is possible to carry out some algebraic operations with Epochs

    # Add 10000 days to a given date
    a = Epoch(1991, 7, 11)
    b = a + 10000
    y, m, d = b.get_date()
    s = str(y) + "/" + str(m) + "/" + str(round(d, 2))
    print_me("1991/7/11 plus 10000 days is", s)

    # Subtract two Epochs to find the number of days between them
    a = Epoch(1986, 2, 9.0)
    b = Epoch(1910, 4, 20.0)
    print_me("The number of days between 1986/2/9 and 1910/4/20 is",
             round(a - b, 2))

    # We can also subtract a given amount of days from an Epoch
    a = Epoch(2003, 12, 31.0)
    b = a - 365.5
    y, m, d = b.get_date()
    s = str(y) + "/" + str(m) + "/" + str(round(d, 2))
    print_me("2003/12/31 minus 365.5 days is", s)

    # Accumulative addition and subtraction of days is also allowed
    a = Epoch(2003, 12, 31.0)
    a += 32.5
    y, m, d = a.get_date()
    s = str(y) + "/" + str(m) + "/" + str(round(d, 2))
    print_me("2003/12/31 plus 32.5 days is", s)

    a = Epoch(2001, 12, 31.0)
    a -= 2 * 365
    y, m, d = a.get_date()
    s = str(y) + "/" + str(m) + "/" + str(round(d, 2))
    print_me("2001/12/31 minus 2*365 days is", s)

    # It is also possible to add days from the right
    a = Epoch(2004, 2, 27.8)
    b = 2.2 + a
    y, m, d = b.get_date()
    s = str(y) + "/" + str(m) + "/" + str(round(d, 2))
    print_me("2004/2/27.8 plus 2.2 days is", s)

    print("")

    # Comparison operadors between epochs are also defined
    a = Epoch(2007, 5, 20.0)
    b = Epoch(2007, 5, 20.000001)
    print_me("2007/5/20.0 == 2007/5/20.000001", a == b)
    print_me("2007/5/20.0 != 2007/5/20.000001", a != b)
    print_me("2007/5/20.0 > 2007/5/20.000001", a > b)
    print_me("2007/5/20.0 <= 2007/5/20.000001", a <= b)

    print("")

    # Compute the time of rise and setting of the Sun in a given day
    e = Epoch(2018, 5, 2)
    print("On May 2nd, 2018, Sun rising/setting times in Munich were (UTC):")
    latitude = Angle(48, 8, 0)
    longitude = Angle(11, 34, 0)
    altitude = 520.0
    rising, setting = e.rise_set(latitude, longitude, altitude)
    y, m, d, h, mi, s = rising.get_full_date()
    print("Rising time: {}:{}".format(h, mi))                           # 3:50
    y, m, d, h, mi, s = setting.get_full_date()
    print("Setting time: {}:{}".format(h, mi))                          # 18:33


if __name__ == "__main__":
    main()
