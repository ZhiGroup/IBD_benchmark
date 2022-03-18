/*
Author: Kecong Tang(Benny)
Group calling module, organize and call core functions.
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace IBD_BM
{
    class GroupCaller:Program
    {
        public static void run(Loader.dataType gtType = Loader.dataType.GroudTruth)
        {
            #region init


            Dictionary<long, List<Loader.IBD_Phy_Start_End>> GT_Holder_All = new Dictionary<long, List<Loader.IBD_Phy_Start_End>>();
            Dictionary<long, List<Loader.IBD_Phy_Start_End>> GT_Holder_gt2cM = new Dictionary<long, List<Loader.IBD_Phy_Start_End>>();

            Dictionary<long, List<Loader.IBD_Phy_Start_End>> RP_Holder_gt2cM = new Dictionary<long, List<Loader.IBD_Phy_Start_End>>();
            Dictionary<long, List<Loader.IBD_Phy_Start_End>> HI_Holder_gt2cM = new Dictionary<long, List<Loader.IBD_Phy_Start_End>>();
            Dictionary<long, List<Loader.IBD_Phy_Start_End>> TP_Holder_gt2cM = new Dictionary<long, List<Loader.IBD_Phy_Start_End>>();
            Dictionary<long, List<Loader.IBD_Phy_Start_End>> IL_Holder_gt2cM = new Dictionary<long, List<Loader.IBD_Phy_Start_End>>();
            Dictionary<long, List<Loader.IBD_Phy_Start_End>> FS_Holder_gt2cM = new Dictionary<long, List<Loader.IBD_Phy_Start_End>>();
            utl.siteDict siteDict = null;

            List<string> RP_ResStr = new List<string>();
            List<string> HI_ResStr = new List<string>();
            List<string> TP_ResStr = new List<string>();
            List<string> IL_ResStr = new List<string>();
            List<string> FS_ResStr = new List<string>();

            for (int i = 0; i < 7; i++)
            {
                RP_ResStr.Add("");
                HI_ResStr.Add("");
                TP_ResStr.Add("");
                IL_ResStr.Add("");
                FS_ResStr.Add("");
            }

            #endregion

            #region load

            utl.GenMapV3 gMap = new utl.GenMapV3(gMap_Path, gMap_PositionCol_Index_ZeroBased, gMap_MapCol_Index_ZeroBased);


            if (parallel_Load == true)
            {
                Parallel.Invoke
                    (
                        () => GT_Holder_gt2cM = Loader.Load_IBD(gtType, GT_Path, 2),
                        () => RP_Holder_gt2cM = Loader.Load_IBD(Loader.dataType.RaPID, RP_Path, 2),
                        () => HI_Holder_gt2cM = Loader.Load_IBD(Loader.dataType.HapIBD, HI_Path, 2),
                        () => TP_Holder_gt2cM = Loader.Load_IBD(Loader.dataType.TPBWT, TP_Path, 2)

                    );
                Parallel.Invoke
                    (

                        () => GT_Holder_All = Loader.Load_IBD(gtType, GT_Path, 0),
                        () => IL_Holder_gt2cM = Loader.Load_IBD(Loader.dataType.iLash, IL_Path, 2),
                        () => FS_Holder_gt2cM = Loader.Load_IBD(Loader.dataType.FastSMC, FS_Path, 2),
                        () => siteDict = new utl.siteDict(vcf_Path)
                    );

            }
            else
            {
                GT_Holder_gt2cM = Loader.Load_IBD(gtType, GT_Path, 2);
                RP_Holder_gt2cM = Loader.Load_IBD(Loader.dataType.RaPID, RP_Path, 2);
                HI_Holder_gt2cM = Loader.Load_IBD(Loader.dataType.HapIBD, HI_Path, 2);
                TP_Holder_gt2cM = Loader.Load_IBD(Loader.dataType.TPBWT, TP_Path, 2);
                GT_Holder_All = Loader.Load_IBD(gtType, GT_Path, 0);
                IL_Holder_gt2cM = Loader.Load_IBD(Loader.dataType.iLash, IL_Path, 2);
                FS_Holder_gt2cM = Loader.Load_IBD(Loader.dataType.FastSMC, FS_Path, 2);
                siteDict = new utl.siteDict(vcf_Path);
            }



            #endregion

            #region compute
            //accuracy

            Parallel.Invoke
                (
                    () => RP_ResStr[0] += (new BinCalculator.Accuracy(GT_Holder_All, RP_Holder_gt2cM, gMap)).Acc_Str,
                    () => HI_ResStr[0] += (new BinCalculator.Accuracy(GT_Holder_All, HI_Holder_gt2cM, gMap)).Acc_Str,
                    () => TP_ResStr[0] += (new BinCalculator.Accuracy(GT_Holder_All, TP_Holder_gt2cM, gMap)).Acc_Str,
                    () => IL_ResStr[0] += (new BinCalculator.Accuracy(GT_Holder_All, IL_Holder_gt2cM, gMap)).Acc_Str,
                    () => FS_ResStr[0] += (new BinCalculator.Accuracy(GT_Holder_All, FS_Holder_gt2cM, gMap)).Acc_Str
                );


            //lenAcc
            Console.WriteLine("Len Accuracy...");

            Parallel.Invoke
                (
                    () => RP_ResStr[1] += (new BinCalculator.LengthAccuracy(GT_Holder_All, RP_Holder_gt2cM, gMap)).Val_Str,
                    () => HI_ResStr[1] += (new BinCalculator.LengthAccuracy(GT_Holder_All, HI_Holder_gt2cM, gMap)).Val_Str,
                    () => TP_ResStr[1] += (new BinCalculator.LengthAccuracy(GT_Holder_All, TP_Holder_gt2cM, gMap)).Val_Str,
                    () => IL_ResStr[1] += (new BinCalculator.LengthAccuracy(GT_Holder_All, IL_Holder_gt2cM, gMap)).Val_Str,
                    () => FS_ResStr[1] += (new BinCalculator.LengthAccuracy(GT_Holder_All, FS_Holder_gt2cM, gMap)).Val_Str
                );


            //LenDis
            Console.WriteLine("Len Disc...");

            Parallel.Invoke
                (
                    () => RP_ResStr[2] += (new BinCalculator.LengthDiscrepancy(GT_Holder_All, RP_Holder_gt2cM, gMap)).Val_Str,
                    () => HI_ResStr[2] += (new BinCalculator.LengthDiscrepancy(GT_Holder_All, HI_Holder_gt2cM, gMap)).Val_Str,
                    () => TP_ResStr[2] += (new BinCalculator.LengthDiscrepancy(GT_Holder_All, TP_Holder_gt2cM, gMap)).Val_Str,
                    () => IL_ResStr[2] += (new BinCalculator.LengthDiscrepancy(GT_Holder_All, IL_Holder_gt2cM, gMap)).Val_Str,
                    () => FS_ResStr[2] += (new BinCalculator.LengthDiscrepancy(GT_Holder_All, FS_Holder_gt2cM, gMap)).Val_Str
                );

            //recall and power
            Console.WriteLine("Single Power...");

            BinCalculator.Power_OneBest RP_Pow = null;
            BinCalculator.Power_OneBest HI_Pow = null;
            BinCalculator.Power_OneBest TP_Pow = null;
            BinCalculator.Power_OneBest IL_Pow = null;
            BinCalculator.Power_OneBest FS_Pow = null;
            Parallel.Invoke
            (
                () => RP_Pow = new BinCalculator.Power_OneBest(GT_Holder_gt2cM, RP_Holder_gt2cM, gMap),
                () => HI_Pow = new BinCalculator.Power_OneBest(GT_Holder_gt2cM, HI_Holder_gt2cM, gMap),
                () => TP_Pow = new BinCalculator.Power_OneBest(GT_Holder_gt2cM, TP_Holder_gt2cM, gMap),
                () => IL_Pow = new BinCalculator.Power_OneBest(GT_Holder_gt2cM, IL_Holder_gt2cM, gMap),
                () => FS_Pow = new BinCalculator.Power_OneBest(GT_Holder_gt2cM, FS_Holder_gt2cM, gMap)
            );
            RP_ResStr[3] += RP_Pow.HitCoverage_Str;
            HI_ResStr[3] += HI_Pow.HitCoverage_Str;
            TP_ResStr[3] += TP_Pow.HitCoverage_Str;
            IL_ResStr[3] += IL_Pow.HitCoverage_Str;
            FS_ResStr[3] += FS_Pow.HitCoverage_Str;

            RP_ResStr[4] += RP_Pow.LengthCoverage_Str;
            HI_ResStr[4] += HI_Pow.LengthCoverage_Str;
            TP_ResStr[4] += TP_Pow.LengthCoverage_Str;
            IL_ResStr[4] += IL_Pow.LengthCoverage_Str;
            FS_ResStr[4] += FS_Pow.LengthCoverage_Str;

            //Multi Power
            //now use full set of reported IBD
            Console.WriteLine("Multi Power...");

            BinCalculator.Power_MultiCoverage RP_mPow = null;
            BinCalculator.Power_MultiCoverage HI_mPow = null;
            BinCalculator.Power_MultiCoverage TP_mPow = null;
            BinCalculator.Power_MultiCoverage IL_mPow = null;
            BinCalculator.Power_MultiCoverage FS_mPow = null;

            Parallel.Invoke
                (
                    () => RP_mPow = new BinCalculator.Power_MultiCoverage(GT_Holder_gt2cM, RP_Holder_gt2cM, gMap, siteDict),
                    () => HI_mPow = new BinCalculator.Power_MultiCoverage(GT_Holder_gt2cM, HI_Holder_gt2cM, gMap, siteDict),
                    () => TP_mPow = new BinCalculator.Power_MultiCoverage(GT_Holder_gt2cM, TP_Holder_gt2cM, gMap, siteDict),
                    () => IL_mPow = new BinCalculator.Power_MultiCoverage(GT_Holder_gt2cM, IL_Holder_gt2cM, gMap, siteDict),
                    () => FS_mPow = new BinCalculator.Power_MultiCoverage(GT_Holder_gt2cM, FS_Holder_gt2cM, gMap, siteDict)
                );
            RP_ResStr[5] += RP_mPow.CoverageHit_Str;
            HI_ResStr[5] += HI_mPow.CoverageHit_Str;
            TP_ResStr[5] += TP_mPow.CoverageHit_Str;
            IL_ResStr[5] += IL_mPow.CoverageHit_Str;
            FS_ResStr[5] += FS_mPow.CoverageHit_Str;

            RP_ResStr[6] += RP_mPow.CoverageCm_Str;
            HI_ResStr[6] += HI_mPow.CoverageCm_Str;
            TP_ResStr[6] += TP_mPow.CoverageCm_Str;
            IL_ResStr[6] += IL_mPow.CoverageCm_Str;
            FS_ResStr[6] += FS_mPow.CoverageCm_Str;

            #endregion


            #region outPut
            StreamWriter sw = new StreamWriter(out_Path);
            for (int i = 0; i < BinCalculator.nBin; i++)
            {
                sw.Write("Bin_" + i + "\t");
            }
            sw.WriteLine();

            List<string> testNames = new List<string>();
            
            testNames.Add("Accuracy");
            testNames.Add("Length Accuracy");
            testNames.Add("Length Discepancy");
            testNames.Add("Recall");
            testNames.Add("Power");
            testNames.Add("Accumulative Recall");
            testNames.Add("Accumulative Power");

            for (int i = 0; i < 7; i++)
            {
                sw.WriteLine(testNames[i]);
                sw.WriteLine("FastSMC\t" + FS_ResStr[i]);
                sw.WriteLine("hapIBD\t" + HI_ResStr[i]);
                sw.WriteLine("iLash\t" + IL_ResStr[i]);
                sw.WriteLine("RaPID\t" + RP_ResStr[i]);
                sw.WriteLine("TPBWT\t" + TP_ResStr[i]);
                sw.WriteLine();
            }

            sw.Close();

            #endregion
        }

    }
}
