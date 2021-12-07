/*
Author: Kecong Tang(Benny)
Calculation module, contains methods for accuracy and power.
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

namespace IBD_Benchmark
{
    class Calculator:Program
    {
        /// <summary>
        /// The number of reported IBD segments that overlap (at least 50%) with a ground true IBD segment.
        /// </summary>
        public class Accuracy
        {
            public long RightCnt = 0;
            public long TotalCnt = 0;
            public double Acc = 0;

            public Accuracy(Dictionary<long, List<Loader.IBD_Phy_Start_End>> GT_Set,
                Dictionary<long, List<Loader.IBD_Phy_Start_End>> TGT_Set, utl.GenMapV3 gMap, double overlapTHD = 0.5)
            {

                Parallel.ForEach(TGT_Set.Keys, (oneKey) =>
                {
                    double aStart, aEnd, bStart, bEnd;
                    foreach (Loader.IBD_Phy_Start_End oneIBD_TGT in TGT_Set[oneKey])
                    {
                        aStart = gMap.getGenLoc(oneIBD_TGT.Start);
                        aEnd = gMap.getGenLoc(oneIBD_TGT.End);

                        Interlocked.Increment(ref TotalCnt);
                        if (GT_Set.ContainsKey(oneKey) == false)
                        {
                            continue;
                        }

                        foreach (Loader.IBD_Phy_Start_End oneIBD_GT in GT_Set[oneKey])
                        {
                            bStart = gMap.getGenLoc(oneIBD_GT.Start);
                            bEnd = gMap.getGenLoc(oneIBD_GT.End);
                            if (utl.rangeOverLap_Percentage(aStart, aEnd, bStart, bEnd) >= overlapTHD)
                            {
                                Interlocked.Increment(ref RightCnt);
                                break;
                            }
                        }

                    }

                });


                Acc = (Convert.ToDouble(RightCnt) / Convert.ToDouble(TotalCnt));
            }
        }

        /// <summary>
        /// Power is the average proportion of true IBD segments by one best reported segment.
        /// </summary>
        public class Power_OneBest
        {
            public double Val = 0;
            public double Cnt = 0;
            public double HitCnt = 0;
            public double HitCoverage = 0;
            double hitTHD = 0.5;
            ConcurrentBag<byte> hitCounter = new ConcurrentBag<byte>();
            ConcurrentBag<double> coverages = new ConcurrentBag<double>();

            public Power_OneBest(Dictionary<long, List<Loader.IBD_Phy_Start_End>> GT_Set,
                Dictionary<long, List<Loader.IBD_Phy_Start_End>> TGT_Set, utl.GenMapV3 gMap, double HitCntCoverageTHD = 0.5)
            {
                hitTHD = HitCntCoverageTHD;
                Parallel.ForEach(GT_Set.Keys, (oneKey) =>
                {
                    double aStart, aEnd, bStart, bEnd;
                    foreach (Loader.IBD_Phy_Start_End oneIBD_GT in GT_Set[oneKey])
                    {
                        aStart = gMap.getGenLoc(oneIBD_GT.Start);
                        aEnd = gMap.getGenLoc(oneIBD_GT.End);


                        if (TGT_Set.ContainsKey(oneKey) == false)
                        {
                            coverages.Add(0);
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
                            coverages.Add(oneCoverages.Max());
                            if (oneCoverages.Max() >= HitCntCoverageTHD)
                            {
                                hitCounter.Add(0);
                            }
                        }
                        else
                        {
                            coverages.Add(0);
                        }
                    }

                });


                Val = coverages.Average();
                Cnt = coverages.Count();
                HitCnt = hitCounter.Count();
                HitCoverage = HitCnt / Cnt;
            }
        }


        /// <summary>
        /// Power is the average proportion of true IBD segments by one best reported segment.
        /// </summary>
        public class Power_MultiCoverage
        {
            public double Coverage_cM = 0;
            public double Coverage_Site = 0;

            ConcurrentBag<double> Lengths = new ConcurrentBag<double>();
            ConcurrentBag<double> LengthCovered = new ConcurrentBag<double>();
            ConcurrentBag<double> nSites = new ConcurrentBag<double>();
            ConcurrentBag<double> nSitesCovered = new ConcurrentBag<double>();

            ConcurrentBag<double> coverages = new ConcurrentBag<double>();

            public Power_MultiCoverage(Dictionary<long, List<Loader.IBD_Phy_Start_End>> GT_Set,
                Dictionary<long, List<Loader.IBD_Phy_Start_End>> TGT_Set, utl.GenMapV3 gMap,
                utl.siteDict siteDict)
            {
                ParallelOptions op = new ParallelOptions();
                op.MaxDegreeOfParallelism = 1;


                Parallel.ForEach(GT_Set.Keys, (oneKey) =>
                {
                    int gtStart, gtEnd, tgtStart, tgtEnd;
                    int Rel_tgtStart, Rel_tgtEnd;//relative index
                    int checkerStart, checkerEnd;

                    BitArray br_Checker;//this checker use relative index

                    foreach (Loader.IBD_Phy_Start_End oneIBD_GT in GT_Set[oneKey])
                    {//for one IBD in GT set
                        gtStart = oneIBD_GT.Start;
                        gtEnd = oneIBD_GT.End;
                        int offSet = siteDict.Get_SiteIndex(gtStart, true);
                        br_Checker = new BitArray(siteDict.Get_SiteIndex(gtEnd, false) - offSet + 1);

                        if (TGT_Set.ContainsKey(oneKey) == false)
                        {//no exist in tgt
                            nSites.Add(br_Checker.Length);
                            Lengths.Add(gMap.getGenDistrance(gtStart, gtEnd));
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


                        Lengths.Add(gMap.getGenDistrance(oneIBD_GT.Start, oneIBD_GT.End));
                        LengthCovered.Add(sumLenCov);
                        //sum up nSite and nSiteCovered
                        int covCnt = 0;
                        foreach (bool one in br_Checker)
                        {
                            if (one == true)
                            {
                                covCnt++;
                            }
                        }


                        nSites.Add(br_Checker.Length);
                        nSitesCovered.Add(covCnt);


                        //end of one GT
                    }

                });
             
                Coverage_cM = LengthCovered.Sum() / Lengths.Sum();
                Coverage_Site = nSitesCovered.Sum() / nSites.Sum();

            }
        }


        /// <summary>
        /// Length Accuracy is the average proportion of reported segments by true IBD segments.
        /// </summary>
        public class LengthAccuracy
        {
            public double Val = 0;
            public double Cnt = 0;


            ConcurrentBag<double> coverages = new ConcurrentBag<double>();

            public LengthAccuracy(Dictionary<long, List<Loader.IBD_Phy_Start_End>> GT_Set,
                Dictionary<long, List<Loader.IBD_Phy_Start_End>> TGT_Set, utl.GenMapV3 gMap)
            {
                Parallel.ForEach(TGT_Set.Keys, (oneKey) =>           
                {
                    //search TGT in GT set
                    double aStart, aEnd, bStart, bEnd;
                    foreach (Loader.IBD_Phy_Start_End oneIBD_TGT in TGT_Set[oneKey])
                    {
                        aStart = gMap.getGenLoc(oneIBD_TGT.Start);
                        aEnd = gMap.getGenLoc(oneIBD_TGT.End);

                        if (GT_Set.ContainsKey(oneKey) == false)
                        {
                            coverages.Add(0);
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
                            coverages.Add(oneCoverages.Max());
                        }
                        else
                        {
                            coverages.Add(0);
                        }
                    }


                });


                Val = coverages.Average();
                Cnt = coverages.Count();

            }
        }


        /// <summary>
        /// Length discrepancy is the mean square root of differences between the lengths of 
        /// ground true segments and (highest overlap) reported segment.
        /// </summary>
        public class LengthDiscrepancy
        {

            public double Val = 0;
            public double Cnt = 0;
            ConcurrentBag<double> squares = new ConcurrentBag<double>();
            public LengthDiscrepancy(Dictionary<long, List<Loader.IBD_Phy_Start_End>> GT_Set,
                Dictionary<long, List<Loader.IBD_Phy_Start_End>> TGT_Set, utl.GenMapV3 gMap)
            {
                Parallel.ForEach(TGT_Set.Keys, (oneKey) =>
                {


                    //search TGT in GT set
                    double aStart, aEnd, bStart, bEnd;
                    foreach (Loader.IBD_Phy_Start_End oneIBD_TGT in TGT_Set[oneKey])
                    {
                        aStart = gMap.getGenLoc(oneIBD_TGT.Start);
                        aEnd = gMap.getGenLoc(oneIBD_TGT.End);


             
                        temCKR.startPhy = oneIBD_TGT.Start;
                        temCKR.endPhy = oneIBD_TGT.End;
                        temCKR.total_Len = aEnd - aStart;
                        temCKR.startGen = aStart;
                        temCKR.endGen = aEnd;
                        temCKR.key = oneKey;

                        if (GT_Set.ContainsKey(oneKey) == false)
                        {
                            squares.Add(Math.Pow(aEnd - aStart, 2));
                            temCKR.coverage = 0;
                            temCKR.dis = Math.Pow(aEnd - aStart, 2);
      
                            continue;
                        }

                        List<double> oneCoverages = new List<double>();

                        //find the best overlap IBD
                        double ovl_Highest = -1;
                        double bStart_Highest = 0, bEnd_Highest = 0;
                        foreach (Loader.IBD_Phy_Start_End oneIBD_GT in GT_Set[oneKey])
                        {
                            bStart = gMap.getGenLoc(oneIBD_GT.Start);
                            bEnd = gMap.getGenLoc(oneIBD_GT.End);
                  
                            double ovl = utl.rangeCovered_Percentage(aStart, aEnd, bStart, bEnd);


                            if (double.IsNaN(ovl) == false && ovl != -1 && ovl > ovl_Highest)
                            {
                                bStart_Highest = bStart;
                                bEnd_Highest = bEnd;
                                ovl_Highest = ovl;
                            }
                        }
                        if (ovl_Highest == -1)
                        {
                            squares.Add(Math.Pow(aEnd - aStart, 2));
                            temCKR.coverage = 0;
                            temCKR.dis = Math.Pow(aEnd - aStart, 2);
                       

                        }
                        else
                        {
                            squares.Add(Math.Pow((aEnd - aStart) - (bEnd_Highest - bStart_Highest), 2));
                            temCKR.coverage = ovl_Highest;
                            temCKR.dis = Math.Pow((aEnd - aStart) - (bEnd_Highest - bStart_Highest), 2);
                            temCKR.matchStart = bStart_Highest;
                            temCKR.mathchEnd = bEnd_Highest;
                        
                        }
                    }
       

                });


                Val = Math.Sqrt(squares.Average());
                Cnt = squares.Count();

            }

        }
    }
}
