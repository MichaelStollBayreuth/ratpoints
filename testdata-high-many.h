/***********************************************************************
 * ratpoints-2.2                                                       *
 *  - A program to find rational points on hyperelliptic curves        *
 * Copyright (C) 2008, 2009, 2022  Michael Stoll                       *
 *                                                                     *
 * This program is free software: you can redistribute it and/or       *
 * modify it under the terms of the GNU General Public License         *
 * as published by the Free Software Foundation, either version 2 of   *
 * the License, or (at your option) any later version.                 *
 *                                                                     *
 * This program is distributed in the hope that it will be useful,     *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of      *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
 * GNU General Public License for more details.                        *
 *                                                                     *
 * You should have received a copy of version 2 of the GNU General     *
 * Public License along with this program.                             *
 * If not, see <http://www.gnu.org/licenses/>.                         *
 ***********************************************************************/

/***********************************************************************
 * testdata-high-many.h                                                *
 *                                                                     *
 * The test curves for "make testhighmany": the thirty curves of       *
 * testdata-many.h that run out of sieving primes.                     *
 *                                                                     *
 * A curve is in this list when sieving_info has to look beyond        *
 * RATPOINTS_DEFAULT_NUM_PRIMES to find enough informative ones, which *
 * is what happens when f is a square modulo every residue for several *
 * of the smallest primes: those say nothing about which numerators    *
 * can occur and are thrown away.  It is the densest curves this       *
 * happens to, so they are the last thirty entries of testdata-many.h, *
 * and they are the ones on which the range of primes is decided by    *
 * the curve rather than by the default.  They span the whole range    *
 * that occurs: between 31 and 36 primes have to be looked at.         *
 *                                                                     *
 * The suite is meant for a large height bound -- see TESTHEIGHT in    *
 * the Makefile, currently 2*10^5 -- where almost all of the time goes *
 * into the sieve.  Thirty curves at that height take about 45         *
 * seconds.  Their answers at 16383 are the last thirty lines of       *
 * testbase-many, so the two lists can be checked against each other.  *
 *                                                                     *
 * Michael Stoll, September 9, 2026                                    *
 ***********************************************************************/

#define NUM_TEST 30

long testdata[NUM_TEST][7] =
{
 {23982047296,-22721787504,-12426175823,13205249358,2410791073,-571054704,29680704},
 {456976,-448032,-255200,208380,61033,-12834,81},
 {265383934425,431854703400,140894275318,-19348834520,-1150028015,-210880560,32901696},
 {4431025,506635,-3907189,-159107,863396,26440,21025},
 {1030276,-2082764,477797,711038,-10667,-24724,1444},
 {2342633056,2577186064,-446912791,-694520420,43527470,10405660,281425},
 {71824,108200,662313,-232358,-382659,59796,49284},
 {57600,-44340,-160751,97015,134217,-101500,18384},
 {27225,718872,3572500,2937044,-160376,-535392,104976},
 {1890625,3407030,-536425,-3162500,136675,529870,75625},
 {166009,172350,-558737,-889176,1083847,25410,441},
 {90440100,-329787420,2108857381,2859496886,399484429,-335264580,396900},
 {1094116,5696856,10995572,-6766946,-859367,462358,37636},
 {1375188363396,1072036985436,134387206573,-6483774878,-145787831,1895880,176400},
 {245238336,-332624544,-685536908,371160148,617289253,100934946,5382369},
 {1942564,2338004,-197633,-568027,47876,15805,1075},
 {-3914316,-9252936,124845003,18011934,-57073455,-10790550,4473225},
 {1737139041,-1007386182,-1162123299,136567300,596158596,-79179328,256},
 {374544,-92520,-981095,392110,823117,10020,900},
 {-96706944,1339600224,10663054192,8322394832,3609958368,3707763494,1490240959},
 {18648633600,-44798214720,-35198667164,81965438100,52229960653,-7147126230,-847565039},
 {15069924,2918256,-31132652,-4234256,18232417,212838,3249},
 {19321,24311876,21561324,-9470570,-406929,178794,21609},
 {479923200,978996480,89706025,-397833150,37690825,12257520,846400},
 {5576200785,2543192226,-7181483345,-1396979652,2440870639,83062434,23261329},
 {9840456000,-8650292160,-35364209936,23301445088,31007000908,-6553516476,438860601},
 {9985600,4412080,-21329999,-2052758,12194777,399776,-37376},
 {50979600,8525160,-269340039,57730464,415483884,-11328660,1587600},
 {102313225,520667840,203415874,-1274482804,-416810599,337611820,116424100},
 {802022400,-4556248320,9004690000,5704997400,-7622962175,1184435670,180906625}
};
