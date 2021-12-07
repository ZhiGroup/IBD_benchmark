/*
Author: Kecong Tang(Benny)
Program entrance, for clear demonstration, only a single example is provied.
Judging each test of IBD results could be running in parallel to save time.
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.IO;

namespace IBD_Benchmark
{
    class Program
    {
        public static utl.siteDict siteDict;
        public static utl.GenMapV3 gMap;

        public static string vcfPath = "";
        public static string gMap_Path = "";
        public static string GT_Path = "";

        public static double IBD_LenTHD = 0;
        public static string log_Path = "";
        static void Main(string[] args)
        {

            GT_Path = "our ground truth path here";
            vcfPath = "input vcf path for IBD software";
            gMap_Path = "your genetic map path here";
            //currently 1st (0 based) column is the physical locations 4th column is the genetic rates
            //you can change column indexs in the utl.GenMapV3 function in singleRun
            singleRun("result of IBD software you want to check", GT_Path, 5, Loader.dataType.RaPID, vcfPath);
        }

        /// <summary>
        /// a single entry to check accuracy and power result for only one tool at a time 
        /// </summary>
        /// <param name="inPath">result from IBD software</param>
        /// <param name="gtPath">ground truth path </param>
        /// <param name="IBD_LenTHD">min IBD length output e.g, 3cM, 5cM, basicly what number given to IBD software</param>
        /// <param name="inType">software name</param>
        /// <param name="vcfPath"></param>
        public static void singleRun(string inPath, string gtPath, double IBD_LenTHD, Loader.dataType inType, string vcfPath)
        {


            Calculator.Accuracy Acc = null;
            Calculator.Power_OneBest Pow = null;
            Calculator.Power_MultiCoverage PowMul = null;
            Calculator.LengthAccuracy LAcc = null;
            Calculator.LengthDiscrepancy LDsc = null;


            gMap = new utl.GenMapV3(gMap_Path, 1, 4);
            siteDict = new utl.siteDict(vcfPath);

            Console.WriteLine("Loading...");
            Dictionary<long, List<Loader.IBD_Phy_Start_End>> resHolder = Loader.Load_IBD(inType, inPath, IBD_LenTHD);
            Dictionary<long, List<Loader.IBD_Phy_Start_End>> GT_Holder = Loader.Load_IBD(Loader.dataType.GroudTruth, gtPath);

            Parallel.Invoke(
                () => resHolder = Loader.Load_IBD(inType, inPath, IBD_LenTHD),
                () => GT_Holder = Loader.Load_IBD(Loader.dataType.GroudTruth, gtPath)
                    );

            StreamWriter sr = new StreamWriter(inPath + ".result");

            Acc = new Calculator.Accuracy(GT_Holder, resHolder, gMap);

            LAcc = new Calculator.LengthAccuracy(GT_Holder, resHolder, gMap);

            LDsc = new Calculator.LengthDiscrepancy(GT_Holder, resHolder, gMap);

            ////for power we need a reload 
            GT_Holder = Loader.Load_IBD(Loader.dataType.GroudTruth, GT_Path, IBD_LenTHD);

            Pow = new Calculator.Power_OneBest(GT_Holder, resHolder, gMap);

            PowMul = new Calculator.Power_MultiCoverage(GT_Holder, resHolder, gMap, siteDict);

            sr.WriteLine(LAcc.Cnt + "\t" + Acc.Acc + "\t" + LAcc.Val + "\t" + LDsc.Val + "\t" + Pow.Val + "\t" + Pow.HitCoverage + "\t" + PowMul.Coverage_cM + "\t" + PowMul.Coverage_Site);


            sr.Close();
            Console.WriteLine(DateTime.Now + " Done.");
        }
    }
}
