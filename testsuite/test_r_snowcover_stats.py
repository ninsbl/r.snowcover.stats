#!/usr/bin/env python3

"""
MODULE:    Test of r.snowcover.stats

AUTHOR(S): Stefan Blumentrath <stbl nve no>

PURPOSE:   Test of r.snowcover.stats

COPYRIGHT: (C) 2022 by NVE, Stefan Blumentrath

This program is free software under the GNU General Public
License (>=v2). Read the file COPYING that comes with GRASS
for details.
"""

import grass.script as gscript

from grass.gunittest.case import TestCase
from grass.gunittest.main import test


class TestBasic(TestCase):
    """The main (and only) test case for the r.snowcover.stats module"""

    # Raster maps be used as inputs (they exist in the NC SPM sample location)
    aoi = "./testsuite/aoi.geojson"
    date_start = "2022-08-22"
    date_end = "2022-08-25"

    @classmethod
    def setUpClass(cls):
        """Ensures expected computational region (and anything else needed)

        These are things needed by all test function but not modified by
        any of them.
        """
        pass
        # We will use specific computational region for our process in case
        # something else is running in parallel with our tests.
        # cls.use_temp_region()
        # Use of of the inputs to set computational region
        # cls.runModule("g.region", raster=cls.test_input_1)

    @classmethod
    def tearDownClass(cls):
        """Remove the temporary region (and anything else we created)"""
        pass
        # cls.del_temp_region()

    def tearDown(self):
        """Remove the output created from the module

        This is executed after each test function run. If we had
        something to set up before each test function run, we would use setUp()
        function.

        Since we remove the raster map after running each test function,
        we can reuse the same name for all the test functions.
        """
        pass
        # self.runModule("g.remove", flags="f", type="raster", name=[self.output])

    def test_output_created(self):
        """Check that the module runs"""
        # run the watershed module
        self.assertModule(
            "r.snowcover.stats",
            aoi=self.aoi,
            date_start=self.date_start,
            date_end=self.date_end,
        )
        # check to see if output is in mapset
        # self.assertRasterExists(self.output, msg="Output was not created")

    def test_missing_parameter(self):
        """Check that the module fails when parameters are missing

        Checks absence of each of the three parameters. Each parameter absence
        is tested separatelly.

        Note that this does not cover all the possible combinations, but it
        tries to simulate most of possible user errors and it should cover
        most of the implementation.
        """
        self.assertModuleFail(
            "r.snowcover.stats",
            msg="The aoi parameter should be required",
        )
        self.assertModuleFail(
            "r.snowcover.stats",
            date_start=self.date_end,
            date_end=self.date_start,
            msg="The start date has to be before end date",
        )

    # def test_output_range(self):
    #     """Check to see if output is within the expected range"""
    #     self.assertModule(
    #         "r.snowcover.stats",
    #         a_input=self.test_input_1,
    #         b_input=self.test_input_2,
    #         output=self.output,
    #     )
    #
    #     min_1, max_1 = get_raster_min_max(self.test_input_1)
    #     min_2, max_2 = get_raster_min_max(self.test_input_2)
    #
    #     reference_min = min_1 + min_2
    #     reference_max = max_1 + max_2
    #
    #     self.assertRasterMinMax(
    #         self.output,
    #         reference_min,
    #         reference_max,
    #         msg="Output exceeds the values computed from inputs",
    #     )


if __name__ == "__main__":
    test()
