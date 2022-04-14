/*
Author: Kecong Tang(Benny)
Core calculation module, contains methods for accuracy and power.
*/
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using System.Collections.Concurrent;
using System.IO;
using System.Collections;

namespace IBD_BM
{
    class BinCalculator:Program
    {

        public static int nBin = (int)((maxBin - minBin) / binLen) + 1;

        static int getBinNum(double cm)
        {

            if (cm >= maxBin)
            {
                return nBin - 1;
            }

            return (int)((cm - minBin) / binLen);

        }

        /// <summary>
        /// The number of reported IBD segments that overlap (at least 50%) with a ground true IBD segment.
        /// </summary>
        public class Accuracy
        {
            public long[] RightCnt = new long[nBin];
            public long[] TotalCnt = new long[nBin];
            public double[] Acc = new double[nBin];

            public string Acc_Str = "";

            public Accuracy(Dictionary<long, List<Loader.IBD_Phy_Start_End>> GT_Set,
                Dictionary<long, List<Loader.IBD_Phy_Start_End>> TGT_Set, utl.GenMapV3 gMap, double overlapTHD = 0.5)
            {

                Parallel.ForEach(TGT_Set.Keys, (oneKey) =>
                //foreach (long oneKey in TGT_Set.Keys)
                {
                    double aStart, aEnd, bStart, bEnd, gLen;
                    int binNum = 0;
                    foreach (Loader.IBD_Phy_Start_End oneIBD_TGT in TGT_Set[oneKey])
                    {
                        aStart = gMap.getGenLoc(oneIBD_TGT.Start);
                        aEnd = gMap.getGenLoc(oneIBD_TGT.End);
                        gLen = gMap.getGenDistrance(oneIBD_TGT.Start, oneIBD_TGT.End);
                        binNum = getBinNum(gLen);
                        Interlocked.Increment(ref TotalCnt[binNum]);
                        if (GT_Set.ContainsKey(oneKey) == false)
                        {
                            continue;
                        }

                        foreach (Loader.IBD_Phy_Start_End oneIBD_GT in GT_Set[oneKey])
                        {
                            bStart = gMap.getGenLoc(oneIBD_GT.Start);
                            bEnd = gMap.getGenLoc(oneIBD_GT.End);
                            if (utl.rangeCovered_Percentage(aStart, aEnd, bStart, bEnd) >= overlapTHD)
                            {
                                Interlocked.Increment(ref RightCnt[binNum]);
                                break;
                            }
                        }

                    }

                });

                for (int i = 0; i < nBin; i++)
                {

                    Acc[i] = (Convert.ToDouble(RightCnt[i]) / Convert.ToDouble(TotalCnt[i]));
                    Acc_Str += Acc[i] + "\t";
                }

            }
        }

        /// <summary>
        /// Power is the average proportion of true IBD segments by one best reported segment.
        /// </summary>
        public class Power_OneBest
        {

            public string HitCoverage_Str = "";
            public string LengthCoverage_Str = "";


            public Power_OneBest(Dictionary<long, List<Loader.IBD_Phy_Start_End>> GT_Set,
                Dictionary<long, List<Loader.IBD_Phy_Start_End>> TGT_Set, utl.GenMapV3 gMap, double HitCntCoverageTHD = 0.5)
            {

                List<ConcurrentBag<byte>> hitCounter = new List<ConcurrentBag<byte>>();
                List<ConcurrentBag<double>> coverages = new List<ConcurrentBag<double>>();
                for (int i = 0; i < nBin; i++)
                {
                    hitCounter.Add(new ConcurrentBag<byte>());
                    coverages.Add(new ConcurrentBag<double>());
                }


                Parallel.ForEach(GT_Set.Keys, (oneKey) =>
                //foreach (long oneKey in GT_Set.Keys)
                {
                    double aStart, aEnd, bStart, bEnd, gLen;
                    int binNum = 0;
                    foreach (Loader.IBD_Phy_Start_End oneIBD_GT in GT_Set[oneKey])
                    {
                        aStart = gMap.getGenLoc(oneIBD_GT.Start);
                        aEnd = gMap.getGenLoc(oneIBD_GT.End);
                        gLen = gMap.getGenDistrance(oneIBD_GT.Start, oneIBD_GT.End);
                        binNum = getBinNum(gLen);

                        if (TGT_Set.ContainsKey(oneKey) == false)
                        {
                            coverages[binNum].Add(0);
                            continue;
                        }
                        //search GT in TGT set
                        List<double> oneCoverages = new List<double>();
                        foreach (Loader.IBD_Phy_Start_End oneIBD_TGT in TGT_Set[oneKey])
                        {
                            bStart = gMap.getGenLoc(oneIBD_TGT.Start);
                            bEnd = gMap.getGenLoc(oneIBD_TGT.End);
                            double cov = utl.rangeCovered_Percentage(aStart, aEnd, bStart, bEnd);
                            if (double.IsNaN(cov) == false)
                            {
                                oneCoverages.Add(cov);

                            }
                        }
                        if (oneCoverages.Count() != 0)
                        {
                            coverages[binNum].Add(oneCoverages.Max());
                            if (oneCoverages.Max() >= HitCntCoverageTHD)
                            {
                                hitCounter[binNum].Add(0);
                            }
                        }
                        else
                        {
                            coverages[binNum].Add(0);
                        }
                    }

                });

                for (int i = 0; i < nBin; i++)
                {
                    if (coverages[i].Count() == 0)
                    {
                        LengthCoverage_Str += "-1\t";
                        HitCoverage_Str += "-1\t";
                    }
                    else
                    {
                        LengthCoverage_Str += coverages[i].Average().ToString() + "\t";
                        HitCoverage_Str += (Convert.ToDouble(hitCounter[i].Count()) / Convert.ToDouble(coverages[i].Count())).ToString() + "\t";
                    }
                }
            }
        }


        /// <summary>
        /// Power is the average proportion of true IBD segments by one best reported segment.
        /// </summary>
        public class Power_MultiCoverage
        {
            public double[] Coverage_cM = new double[nBin];
            public double[] Coverage_Site = new double[nBin];
            public double[] Coverage_Hit = new double[nBin];
            public string CoverageCm_Str = "";
            public string CoverageSite_Str = "";
            public string CoverageHit_Str = "";

            public Power_MultiCoverage(Dictionary<long, List<Loader.IBD_Phy_Start_End>> GT_Set,
                Dictionary<long, List<Loader.IBD_Phy_Start_End>> TGT_Set, utl.GenMapV3 gMap,
                utl.siteDict siteDict)
            {

                List<ConcurrentBag<double>> Lengths = new List<ConcurrentBag<double>>();
                List<ConcurrentBag<double>> LengthCovered = new List<ConcurrentBag<double>>();
                List<ConcurrentBag<double>> nSites = new List<ConcurrentBag<double>>();
                List<ConcurrentBag<double>> nSitesCovered = new List<ConcurrentBag<double>>();
                List<ConcurrentBag<double>> coverages = new List<ConcurrentBag<double>>();
                List<ConcurrentBag<byte>> hitCounter = new List<ConcurrentBag<byte>>();
                for (int i = 0; i < nBin; i++)
                {
                    Lengths.Add(new ConcurrentBag<double>());
                    LengthCovered.Add(new ConcurrentBag<double>());
                    nSites.Add(new ConcurrentBag<double>());
                    nSitesCovered.Add(new ConcurrentBag<double>());
                    coverages.Add(new ConcurrentBag<double>());
                    hitCounter.Add(new ConcurrentBag<byte>());
                }



                Parallel.ForEach(GT_Set.Keys, (oneKey) =>
                //foreach (long oneKey in GT_Set.Keys)
                {
                    int gtStart, gtEnd, tgtStart, tgtEnd;
                    int Rel_tgtStart, Rel_tgtEnd;//relative index
                    int checkerStart, checkerEnd;

                    double gLen = 0;
                    int binNum = 0;

                    BitArray br_Checker;//this checker use relative index

                    foreach (Loader.IBD_Phy_Start_End oneIBD_GT in GT_Set[oneKey])
                    {//for one IBD in GT set
                        gtStart = oneIBD_GT.Start;
                        gtEnd = oneIBD_GT.End;
                        int offSet = siteDict.Get_SiteIndex(gtStart, true);
                        br_Checker = new BitArray(siteDict.Get_SiteIndex(gtEnd, false) - offSet + 1);

                        gLen = gMap.getGenDistrance(oneIBD_GT.Start, oneIBD_GT.End);
                        binNum = getBinNum(gLen);
                        Lengths[binNum].Add(gMap.getGenDistrance(gtStart, gtEnd));

                        if (TGT_Set.ContainsKey(oneKey) == false)
                        {//no exist in tgt
                            nSites[binNum].Add(br_Checker.Length);

                            continue;
                        }
                        //search GT in TGT set, update checker                       
                        foreach (Loader.IBD_Phy_Start_End oneIBD_TGT in TGT_Set[oneKey])
                        {
                            tgtStart = oneIBD_TGT.Start;
                            tgtEnd = oneIBD_TGT.End;
                            //update checker


                            Rel_tgtStart = siteDict.Get_SiteIndex(tgtStart, true) - offSet;
                            Rel_tgtEnd = siteDict.Get_SiteIndex(tgtEnd, false) - offSet;
                            checkerStart = Math.Max(Rel_tgtStart, 0);
                            checkerEnd = Math.Min(Rel_tgtEnd, br_Checker.Length - 1);




                            for (int i = checkerStart; i <= checkerEnd; i++)
                            {
                                br_Checker[i] = true;
                            }
                        }

                        //sum up length and length coverage
                        double sumLenCov = 0;
                        int covStart = 0, covEnd = 0;
                        for (int i = 0; i < br_Checker.Length; i++)
                        {
                            if (br_Checker[i] == true)
                            {//starting
                                if (i == 0 || br_Checker[i - 1] == false)
                                {
                                    covStart = i;
                                }
                            }
                            else
                            {//ending
                                if (i != 0 && br_Checker[i - 1] == true)
                                {
                                    covEnd = i - 1;
                                    //update
                                    sumLenCov += gMap.getGenDistrance(siteDict.Get_Phy(covStart + offSet), siteDict.Get_Phy(covEnd + offSet));
                                }

                            }
                        }
                        if (br_Checker[br_Checker.Length - 1] == true)
                        {
                            covEnd = br_Checker.Length - 1;
                            sumLenCov += gMap.getGenDistrance(siteDict.Get_Phy(covStart + offSet), siteDict.Get_Phy(covEnd + offSet));
                        }


                        LengthCovered[binNum].Add(sumLenCov);



                        //update hit counter
                        if (sumLenCov >= gLen / 2)
                        {
                            hitCounter[binNum].Add(0);
                        }


                        //sum up nSite and nSiteCovered
                        int covCnt = 0;
                        foreach (bool one in br_Checker)
                        {
                            if (one == true)
                            {
                                covCnt++;
                            }
                        }


                        nSites[binNum].Add(br_Checker.Length);
                        nSitesCovered[binNum].Add(covCnt);


                        //end of one GT
                    }

                });
                //}
                for (int i = 0; i < nBin; i++)
                {
                    if (LengthCovered[i].Count() == 0 || Lengths[i].Count() == 0)
                    {
                        Coverage_cM[i] = -1;
                    }
                    else
                    {
                        Coverage_cM[i] = LengthCovered[i].Sum() / Lengths[i].Sum();
                    }

                    if (nSitesCovered[i].Count() == 0 || nSites[i].Count() == 0)
                    {
                        Coverage_Site[i] = -1;
                    }
                    else
                    {
                        Coverage_Site[i] = nSitesCovered[i].Sum() / nSites[i].Sum();
                    }

                    if (Lengths[i].Count() == 0)
                    {
                        Coverage_Hit[i] = -1;
                    }
                    else
                    {
                        Coverage_Hit[i] = Convert.ToDouble(hitCounter[i].Count()) / Convert.ToDouble(Lengths[i].Count());
                    }


                    CoverageCm_Str += Coverage_cM[i] + "\t";
                    CoverageSite_Str += Coverage_Site[i] + "\t";
                    CoverageHit_Str += Coverage_Hit[i] + "\t";
                }


            }
        }


        /// <summary>
        /// Length Accuracy is the average proportion of reported segments by true IBD segments.
        /// </summary>
        public class LengthAccuracy
        {
            public double[] Val = new double[nBin];
            public double[] Cnt = new double[nBin];

            public string Val_Str = "";
            public string Cnt_Str = "";

            public LengthAccuracy(Dictionary<long, List<Loader.IBD_Phy_Start_End>> GT_Set,
                Dictionary<long, List<Loader.IBD_Phy_Start_End>> TGT_Set, utl.GenMapV3 gMap)
            {

                List<ConcurrentBag<double>> coverages = new List<ConcurrentBag<double>>();
                for (int i = 0; i < nBin; i++)
                {
                    coverages.Add(new ConcurrentBag<double>());
                }


                Parallel.ForEach(TGT_Set.Keys, (oneKey) =>
                //foreach (long oneKey in TGT_Set.Keys)
                {
                    //search TGT in GT set
                    double aStart, aEnd, bStart, bEnd;
                    double gLen = 0;
                    int binNum = 0;
                    foreach (Loader.IBD_Phy_Start_End oneIBD_TGT in TGT_Set[oneKey])
                    {
                        aStart = gMap.getGenLoc(oneIBD_TGT.Start);
                        aEnd = gMap.getGenLoc(oneIBD_TGT.End);
                        gLen = gMap.getGenDistrance(oneIBD_TGT.Start, oneIBD_TGT.End);
                        binNum = getBinNum(gLen);
                        if (GT_Set.ContainsKey(oneKey) == false)
                        {
                            coverages[binNum].Add(0);
                            continue;
                        }

                        List<double> oneCoverages = new List<double>();
                        foreach (Loader.IBD_Phy_Start_End oneIBD_GT in GT_Set[oneKey])
                        {
                            bStart = gMap.getGenLoc(oneIBD_GT.Start);
                            bEnd = gMap.getGenLoc(oneIBD_GT.End);
                            double cov = utl.rangeCovered_Percentage(aStart, aEnd, bStart, bEnd);
                            if (double.IsNaN(cov) == false)
                            {
                                oneCoverages.Add(cov);
                            }
                        }

                        if (oneCoverages.Count() != 0)
                        {
                            coverages[binNum].Add(oneCoverages.Max());
                        }
                        else
                        {
                            coverages[binNum].Add(0);
                        }
                    }

                    //}

                });

                for (int i = 0; i < nBin; i++)
                {
                    if (coverages[i].Count() == 0)
                    {
                        Val[i] = -1;
                        Cnt[i] = -1;
                    }
                    else
                    {
                        Val[i] = coverages[i].Average();
                        Cnt[i] = coverages[i].Count();
                    }

                    Val_Str += Val[i] + "\t";
                    Cnt_Str += Cnt[i] + "\t";

                }

            }
        }

        /// <summary>
        /// Length discrepancy is the root mean square of differences between the lengths of 
        /// ground true segments and (highest overlap) reported segment.
        /// Note, the max discrepancy should be length of reported IBD itself, will not look for higher coverage but high discrepancy.
        /// E.g:
        /// old logic:  a 10cM GT covers   2cM reported 100%; a 1.5 cM covers this 2cM reported 75%. then this 10cM will be chosen, which is bad.
        /// new logic: set max discrepancy = reported itself, find all overlap GT, then take the min discrepancy
        /// </summary>
        public class LengthDiscrepancy
        {


            public double[] Val = new double[nBin];
            public double[] Cnt = new double[nBin];

            public string Val_Str = "";
            public string Cnt_Str = "";


            public LengthDiscrepancy(Dictionary<long, List<Loader.IBD_Phy_Start_End>> GT_Set,
                Dictionary<long, List<Loader.IBD_Phy_Start_End>> TGT_Set, utl.GenMapV3 gMap)
            {

                List<ConcurrentBag<double>> squares = new List<ConcurrentBag<double>>();
                for (int i = 0; i < nBin; i++)
                {
                    squares.Add(new ConcurrentBag<double>());
                }
                Parallel.ForEach(TGT_Set.Keys, (oneKey) =>
                //foreach (long oneKey in TGT_Set.Keys)
                {
                    //search TGT in GT set
                    double aStart, aEnd, bStart, bEnd;
                    double gLen = 0;
                    int binNum = 0;

                    foreach (Loader.IBD_Phy_Start_End oneIBD_TGT in TGT_Set[oneKey])
                    {
                        aStart = gMap.getGenLoc(oneIBD_TGT.Start);
                        aEnd = gMap.getGenLoc(oneIBD_TGT.End);
                        gLen = gMap.getGenDistrance(oneIBD_TGT.Start, oneIBD_TGT.End);
                        binNum = getBinNum(gLen);

                        double squaresDiff = Math.Pow(gLen, 2);
                        double temSquaresDiff = double.MaxValue;
                        if (GT_Set.ContainsKey(oneKey) == false)
                        {
                            squares[binNum].Add(squaresDiff);
                            continue;
                        }

                        //find the min disc 

                        foreach (Loader.IBD_Phy_Start_End oneIBD_GT in GT_Set[oneKey])
                        {
                            bStart = gMap.getGenLoc(oneIBD_GT.Start);
                            bEnd = gMap.getGenLoc(oneIBD_GT.End);

                            double ovl = utl.rangeCovered_Percentage(aStart, aEnd, bStart, bEnd);
                            if (double.IsNaN(ovl) == false && ovl != 0)
                            {
                                temSquaresDiff = Math.Pow((aEnd - aStart) - (bEnd - bStart), 2);
                                if (temSquaresDiff < squaresDiff)
                                {
                                    squaresDiff = temSquaresDiff;
                                }
                            }
                        }

                        squares[binNum].Add(squaresDiff);

                    }
                    //}

                });

                for (int i = 0; i < nBin; i++)
                {
                    if (squares[i].Count() == 0)
                    {
                        Val[i] = -1;
                        Cnt[i] = -1;
                    }
                    else
                    {
                        Val[i] = Math.Sqrt(squares[i].Average());
                        Cnt[i] = squares[i].Count();
                    }

                    Val_Str += Val[i] + "\t";
                    Cnt_Str += Cnt[i] + "\t";
                }

            }

        }
    }
}
