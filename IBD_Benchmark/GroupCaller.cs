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
            Dictionary<long, List<Loader.IBD_Phy_Start_End>> GT_Holder_TGT_Len = new Dictionary<long, List<Loader.IBD_Phy_Start_End>>();

            Dictionary<long, List<Loader.IBD_Phy_Start_End>> Reported_Holder_TGT_Len = new Dictionary<long, List<Loader.IBD_Phy_Start_End>>();

            utl.siteDict siteDict = null;

            List<string> Reported_ResStr = new List<string>();


            for (int i = 0; i < 7; i++)
            {
                Reported_ResStr.Add("");

            }

            #endregion

            #region load

            utl.GenMapV3 gMap = new utl.GenMapV3(gMap_Path, gMap_PositionCol_Index_ZeroBased, gMap_MapCol_Index_ZeroBased);
            GT_Holder_All = Loader.Load_IBD(gtType, GT_Path, 0);
            GT_Holder_TGT_Len = Loader.Load_IBD(gtType, GT_Path, minBin);
            Reported_Holder_TGT_Len = Loader.Load_IBD(tool_Type, reported_Path, minBin);
            siteDict = new utl.siteDict(vcf_Path);


            #endregion

            #region compute
            //accuracy
            Console.WriteLine("Accuracy...");
            Reported_ResStr[0] += (new BinCalculator.Accuracy(GT_Holder_All, Reported_Holder_TGT_Len, gMap)).Acc_Str;

            //lenAcc
            Console.WriteLine("Len Accuracy...");
            Reported_ResStr[1] += (new BinCalculator.LengthAccuracy(GT_Holder_All, Reported_Holder_TGT_Len, gMap)).Val_Str;

            //LenDis
            Console.WriteLine("Len Disc...");

            Reported_ResStr[2] += (new BinCalculator.LengthDiscrepancy(GT_Holder_All, Reported_Holder_TGT_Len, gMap)).Val_Str;

            //recall and power
            Console.WriteLine("Single Power...");

            BinCalculator.Power_OneBest reported_Pow = new BinCalculator.Power_OneBest(GT_Holder_TGT_Len, Reported_Holder_TGT_Len, gMap);

            Reported_ResStr[3] += reported_Pow.HitCoverage_Str;

            Reported_ResStr[4] += reported_Pow.LengthCoverage_Str;

            //Multi Power
            //now use full set of reported IBD
            Console.WriteLine("Multi Power...");

            BinCalculator.Power_MultiCoverage reported_mPow = new BinCalculator.Power_MultiCoverage(GT_Holder_TGT_Len, Reported_Holder_TGT_Len, gMap, siteDict);

            Reported_ResStr[5] += reported_mPow.CoverageHit_Str;

            Reported_ResStr[6] += reported_mPow.CoverageCm_Str;


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
                sw.WriteLine(tool_Name+"\t"+Reported_ResStr[i]);
                sw.WriteLine();
            }

            sw.Close();

            #endregion
        }

    }
}
