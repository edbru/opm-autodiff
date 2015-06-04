/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <config.h>

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE WellStateResizeTest

#include <opm/core/wells.h>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/grid/GridManager.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/core/simulator/BlackoilState.hpp>

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <memory>

using namespace Opm;
/*
Example of COORD and ZCORN keywords found in opm-parser/opm/parser/eclipse/EclipseState/Grid/tests/EclipseGridTests.cpp.
input2 (with DXV, DYV, DZV and TOPS keywords) found in opm-core/tests/test_EclipseWriteRFTHandler.cpp.
If problem with the input data deck, the file opm-core/tests/wells_manager_data.data seems to have data that works
(gives us wells in c_wells() in the wellsManager).
*/

void testSizes(EclipseStateConstPtr eclipseState, WellState& wsResize, Opm::GridManager &gridManager, int reportStep){
    wsResize.resize(*eclipseState, reportStep);

    int numBhpResize = wsResize.bhp().size();
    int numTemperatureResize = wsResize.temperature().size();
    int numWellRatesResize = wsResize.wellRates().size();
    int numPerfRatesResize = wsResize.perfRates().size();
    int numPerfPressResize = wsResize.perfPress().size();

    Opm::WellsManager wellsManager(eclipseState, reportStep , *gridManager.c_grid(), NULL);
    const Wells* wells = wellsManager.c_wells();
    std::shared_ptr<Opm::GridManager> ourFineGridManagerPtr(new Opm::GridManager(eclipseState->getEclipseGrid()));
    std::shared_ptr<Opm::BlackoilState> blackoilState(new Opm::BlackoilState);
    blackoilState->init(*ourFineGridManagerPtr->c_grid(), 3);

    BlackoilState* st = blackoilState.get();
    WellState wsInit;
    wsInit.init(wells, *st);
    int numBhpInit = wsInit.bhp().size();
    int numTemperatureInit = wsInit.temperature().size();
    int numWellRatesInit = wsInit.wellRates().size();
    int numPerfRatesInit = wsInit.perfRates().size();
    int numPerfPressInit = wsInit.perfPress().size();

    BOOST_CHECK((numBhpResize == numBhpInit) &&
                (numTemperatureResize == numTemperatureInit) &&
                (numWellRatesResize == numWellRatesInit) &&
                (numPerfRatesResize == numPerfRatesInit) &&
                (numPerfPressResize == numPerfPressInit));
}

BOOST_AUTO_TEST_CASE(resizeWellState) {
    Opm::Parser parser;
/*
    std::string input =
            "START             -- 0 \n"
            "19 JUN 2007 / \n"

            "RUNSPEC\n"

            "OIL\n"
            "GAS\n"
            "WATER\n"
            "DIMENS\n"
            " 10 10 10 /\n"

            "GRID\n"

            "DXV\n"
            "10*1000.0 /\n"

            "DYV\n"
            "10*1000.0 /\n"

            "DZV\n"
            "10.0 20.0 30.0 10.0 5.0 /\n"

            "TOPS\n"
            "100*10 /\n"

            "COORD\n"
            "  726*1 / \n"
            "ZCORN \n"
            "  8000*1 / \n"
            "SCHEDULE\n"

            "DATES             -- 1\n"
            " 10  OKT 2008 / \n"
            "/\n"
            "WELSPECS\n"
            "    'OP_1'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
            "    'OP_2'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
            "/\n"
            "COMPDAT\n"
            " 'OP_1'  9  9   1   1 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            " 'OP_2'  9  9   2   2 'OPEN' 1*   46.825   0.311  4332.346 1*  1*  'X'  22.123 / \n"
            " 'OP_1'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            " 'OP_1'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            "/\n"

            "DATES             -- 2\n"
            " 20  JAN 2010 / \n"
            "/\n"
            "WELSPECS\n"
            "    'OP_3'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
            "/\n"

            "DATES             -- 3\n"
            " 15  JUN 2013 / \n"
            "/\n"
            "COMPDAT\n"
            " 'OP_2'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            " 'OP_1'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            "/\n"

            "DATES             -- 4\n"
            " 22  APR 2014 / \n"
            "/\n"
            "WELSPECS\n"
            "    'OP_4'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
            "/\n"
            "COMPDAT\n"
            " 'OP_2'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            " 'OP_1'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            " 'OP_2'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            " 'OP_3'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            " 'OP_4'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            " 'OP_3'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            " 'OP_2'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            "/\n";
*/

/*
    std::string input1 =
                "RUNSPEC\n"
                "OIL\n"
                "GAS\n"
                "WATER\n"
                "DIMENS\n"
                " 10 10 10 /\n"
                "GRID\n"
                "DXV\n"
                "10*0.25 /\n"
                "DYV\n"
                "10*0.25 /\n"
                "DZV\n"
                "10*0.25 /\n"
                "TOPS\n"
                "100*0.25 /\n"
                "\n"
                 "START             -- 0 \n"
                "1 NOV 1979 / \n"
                "SCHEDULE\n"
                "DATES             -- 1\n"
                " 1 DES 1979/ \n"
                "/\n"
                "WELSPECS\n"
                "    'OP_1'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
                "    'OP_2'       'OP'   4   4 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
                "/\n"
                "COMPDAT\n"
                " 'OP_1'  9  9   1   1 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
                " 'OP_1'  9  9   2   2 'OPEN' 1*   46.825   0.311  4332.346 1*  1*  'X'  22.123 / \n"
                " 'OP_1'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
                " 'OP_2'  4  4   4  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
                "/\n"
                "DATES             -- 2\n"
                " 10  OKT 2008 / \n"
                "/\n"
                "WRFT \n"
                "/ \n"
                "WELOPEN\n"
                " 'OP_1' OPEN / \n"
                " 'OP_2' OPEN / \n"
                "/\n"
                "DATES             -- 3\n"
                " 10  NOV 2008 / \n"
                "/\n";
*/


    std::string input2 =
            "RUNSPEC\n"
            "OIL\n"
            "GAS\n"
            "WATER\n"
            "DIMENS\n"
            " 10 10 10 /\n"

            "GRID\n"
            "DXV\n"
            "10*0.25 /\n"
            "DYV\n"
            "10*0.25 /\n"
            "DZV\n"
            "10*0.25 /\n"
            "TOPS\n"
            "100*0.25 /\n"
            "\n"

            "START             -- 0 \n"
            "1 NOV 1979 / \n"

            "SCHEDULE\n"
            "DATES             -- 1\n"
            " 10  OKT 2008 / \n"
            "/\n"
            "WELSPECS\n"
            "    'OP_1'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
            "    'OP_2'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
            "/\n"
            "COMPDAT\n"
            " 'OP_1'  9  9   1   1 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            " 'OP_2'  9  9   2   2 'OPEN' 1*   46.825   0.311  4332.346 1*  1*  'X'  22.123 / \n"
            " 'OP_1'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            "/\n"
            "WCONPROD\n"
                "'OP_1' 'OPEN' 'ORAT' 20000  4* 1000 /\n"
            "/\n"
            "WCONINJE\n"
                "'OP_2' 'GAS' 'OPEN' 'RATE' 100 200 400 /\n"
            "/\n"

            "DATES             -- 2\n"
            " 20  JAN 2011 / \n"
            "/\n"
            "WELSPECS\n"
            "    'OP_3'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
            "/\n"
            "COMPDAT\n"
            " 'OP_3'  9  9   1   1 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            "/\n"
            "WCONPROD\n"
                "'OP_3' 'OPEN' 'ORAT' 20000  4* 1000 /\n"
            "/\n"

            "DATES             -- 3\n"
            " 15  JUN 2013 / \n"
            "/\n"
            "COMPDAT\n"
            " 'OP_2'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            " 'OP_1'  9  9   7  7 'SHUT' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            "/\n"

            "DATES             -- 4\n"
            " 22  APR 2014 / \n"
            "/\n"
            "WELSPECS\n"
            "    'OP_4'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
            "/\n"
            "COMPDAT\n"
            " 'OP_2'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            " 'OP_4'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            " 'OP_3'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            "/\n"
            "WCONPROD\n"
                "'OP_4' 'OPEN' 'ORAT' 20000  4* 1000 /\n"
            "/\n"

            "DATES             -- 5\n"
            " 30  AUG 2014 / \n"
            "/\n"
            "WELSPECS\n"
            "    'OP_5'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
            "/\n"
            "COMPDAT\n"
            " 'OP_5'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            "/\n"
            "WCONPROD\n"
                "'OP_5' 'OPEN' 'ORAT' 20000  4* 1000 /\n"
            "/\n"

            "DATES             -- 6\n"
            " 15  SEP 2014 / \n"
            "/\n"
            "WCONPROD\n"
                "'OP_3' 'SHUT' 'ORAT' 20000  4* 1000 /\n"
            "/\n"

            "DATES             -- 7\n"
            " 9  OCT 2014 / \n"
            "/\n"
            "WELSPECS\n"
            "    'OP_6'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
            "/\n"
            "COMPDAT\n"
            " 'OP_6'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            "/\n"
            "WCONPROD\n"
                "'OP_6' 'OPEN' 'ORAT' 20000  4* 1000 /\n"
            "/\n";

/*
    std::string input3 =
    "RUNSPEC\n"

    "OIL\n"
    "GAS\n"
    "WATER\n"

    "DIMENS\n"
    " 10 10  5  /\n"

    "GRID\n"

    "DXV\n"
    "10*1000.0 /\n"

    "DYV\n"
    "10*1000.0 /\n"

    "DZV\n"
    "10.0 20.0 30.0 10.0 5.0 /\n"

    "TOPS\n"
    " 100*10 /\n"

    "SCHEDULE\n"

    "WELSPECS\n"
        "'INJ1' 'G'    1  1    8335 'GAS'  /\n"
        "'PROD1' 'G'   10 10    8400 'OIL'  /\n"
    "/\n"

    "COMPDAT\n"
        "'INJ1'   1  1 1  1 'OPEN' 1   10.6092   0.5  /\n"
        "'PROD1'  10 3 3  3 'OPEN' 0   10.6092   0.5  /\n"
    "/\n"

    "WCONPROD\n"
         "'PROD1' 'OPEN' 'ORAT' 20000  4* 1000 /\n"
    "/\n"

    "WCONINJE\n"
         "'INJ1' 'GAS' 'OPEN' 'RATE' 100 200 400 /\n"
    "/\n"


    "DATES   -- Step1\n"
     "1 'FEB' 2000 /\n"
    "/\n"

    "WCONPROD\n"
       "'PROD1' 'OPEN' 'RESV' 999  3* 123 100 /\n"
    "/\n"

    "WCONINJE\n"
       "'INJ1' 'WATER' 'OPEN' 'RESV' 10 20 40 /\n"
    "/\n"


    "DATES  -- Step2\n"
     "1 'MAR' 2000 /\n"
    "/\n"


    "WCONPROD\n"
      "'PROD1'  'SHUT' /\n"
    "/\n"


    "DATES  -- Step3\n"
     "1 'APR' 2000 /\n"
    "/\n"

    "WELSPECS\n"
       "'NEW'  'G'   2   2  1*       'OIL'  2*      'STOP'  4* /\n"
    "/\n"

    "COMPDAT\n"
        "'NEW'   2  2 2  2 'OPEN' 1   10.6092   0.5  /\n"
    "/\n"


    "WCONHIST\n"
       "'NEW'      'OPEN'      'ORAT'      0.000      0.000      0.000  5* /\n"
    "/\n"

    "END\n";
*/

    double i1BhpSet = 1.25;
    double i2BhpSet = 2.35;
    double i3BhpSet = 3.45;
    double i1TemperatureSet = 4.55;
    double i2TemperatureSet = 5.65;
    double i3TemperatureSet = 6.75;
    double i1WellRatesSet = 7.25;
    double i2WellRatesSet = 8.35;
    double i3WellRatesSet = 9.05;
    double i1PerfRatesSet = 11.25;
    double i2PerfRatesSet = 12.25;
    double i3PerfRatesSet = 13.25;
    double i1PerfPressSet = 14.25;
    double i2PerfPressSet = 15.25;
    double i3PerfPressSet = 16.25;

    DeckConstPtr deck = parser.parseString(input2);
    EclipseStateConstPtr eclipseState = std::make_shared<const EclipseState>(deck);
    Opm::GridManager gridManager(deck);

    WellState wsResize;
    testSizes(eclipseState, wsResize, gridManager, 1);

    wsResize.setBhpValue(i1BhpSet, "OP_1");
    wsResize.setBhpValue(i2BhpSet, 1);
    wsResize.setTemperatureValue(i1TemperatureSet, "OP_1");
    wsResize.setTemperatureValue(i2TemperatureSet, 1);
    wsResize.setWellRatesValue(i1WellRatesSet, "OP_1", 1);
    wsResize.setWellRatesValue(i2WellRatesSet, 5);
    wsResize.setPerfRatesValue(i1PerfRatesSet, 0);
    wsResize.setPerfRatesValue(i2PerfRatesSet, 1);
    wsResize.setPerfPressValue(i1PerfPressSet, 0);
    wsResize.setPerfPressValue(i2PerfPressSet, 1);
    wsResize.setPerfRatesValue(i3PerfRatesSet, 7);
    wsResize.setPerfPressValue(i3PerfPressSet, 7);

    testSizes(eclipseState, wsResize, gridManager, 2);
    testSizes(eclipseState, wsResize, gridManager, 3);
    testSizes(eclipseState, wsResize, gridManager, 4);

    wsResize.setBhpValue(i3BhpSet, 3);
    wsResize.setTemperatureValue(i3TemperatureSet, 3);
    wsResize.setWellRatesValue(i3WellRatesSet, 9);

    testSizes(eclipseState, wsResize, gridManager, 5);
    testSizes(eclipseState, wsResize, gridManager, 6);
    testSizes(eclipseState, wsResize, gridManager, 7);

    double i1BhpGet = wsResize.getBhpValue(0);
    double i2BhpGet = wsResize.getBhpValue("OP_2");
    double i3BhpGet = wsResize.getBhpValue(2);
    double i1TemperatureGet = wsResize.getTemperatureValue(0);
    double i2TemperatureGet = wsResize.getTemperatureValue("OP_2");
    double i3TemperatureGet = wsResize.getTemperatureValue(2);
    double i1WellRatesGet = wsResize.getWellRatesValue(0);
    double i2WellRatesGet = wsResize.getWellRatesValue("OP_2", 3);
    double i3WellRatesGet = wsResize.getWellRatesValue(6);
    double i1PerfRatesGet = wsResize.getPerfRatesValue(0);
    double i2PerfRatesGet = wsResize.getPerfRatesValue(1);
    double i3PerfRatesGet = wsResize.getPerfRatesValue(6);
    double i1PerfPressGet = wsResize.getPerfPressValue(0);
    double i2PerfPressGet = wsResize.getPerfPressValue(1);
    double i3PerfPressGet = wsResize.getPerfPressValue(6);

    BOOST_CHECK_EQUAL(i1BhpGet, i1BhpSet);
    BOOST_CHECK_EQUAL(i2BhpGet, i2BhpSet);
    BOOST_CHECK_EQUAL(i3BhpGet, i3BhpSet);
    BOOST_CHECK_EQUAL(i1TemperatureGet, i1TemperatureSet);
    BOOST_CHECK_EQUAL(i2TemperatureGet, i2TemperatureSet);
    BOOST_CHECK_EQUAL(i3TemperatureGet, i3TemperatureSet);
    BOOST_CHECK_EQUAL(i1WellRatesGet, i1WellRatesSet);
    BOOST_CHECK_EQUAL(i2WellRatesGet, i2WellRatesSet);
    BOOST_CHECK_EQUAL(i3WellRatesGet, i3WellRatesSet);
    BOOST_CHECK_EQUAL(i1PerfRatesGet, i1PerfRatesSet);
    BOOST_CHECK_EQUAL(i2PerfRatesGet, i2PerfRatesSet);
    BOOST_CHECK_EQUAL(i3PerfRatesGet, i3PerfRatesSet);
    BOOST_CHECK_EQUAL(i1PerfPressGet, i1PerfPressSet);
    BOOST_CHECK_EQUAL(i2PerfPressGet, i2PerfPressSet);
    BOOST_CHECK_EQUAL(i3PerfPressGet, i3PerfPressSet);
}
