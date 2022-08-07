using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using System.IO;

namespace IBD_BM
{
    class siteCoverage
    {

        /// <summary>
        /// generate a IBD coverage data by site
        /// each site has a counter
        /// this module requires: a vcf file to load site dictionary 
        /// 
        /// Usage:
        /// 1. Assign vcf_Path, IBD_result_Path ( results of IBD detection tools), and output path
        /// 2. Choose IBD data type: in Loader.Load_IBD module.
        /// 3. Call singleRun() in program.cs/ Main to run.
        /// </summary>
        public static void singleRun()
        {
            string vcf_Path = "E:\\s15in\\mix.arr.e0.001.vcf";
            string IBD_result_Path = "E:\\s15in\\mix.gt.arr.txt";
            string outPath="E:\\tem\\dummy.out.txt";

            Dictionary<long, List<Loader.IBD_Phy_Start_End>> IBDs = Loader.Load_IBD(Loader.dataType.GroudTruth, IBD_result_Path);
            utl.siteDict siteD = new utl.siteDict(vcf_Path);

            runSiteCoverage(IBDs, siteD, outPath);


        }

        static List<int> runSiteCoverage(Dictionary<long, List<Loader.IBD_Phy_Start_End>> IBDs, utl.siteDict sitDict, string outPath)
        {

            int[] cnts = new int[sitDict.Count()];
            ParallelOptions op = new ParallelOptions();
            op.MaxDegreeOfParallelism = 10;

            Parallel.ForEach(IBDs.Keys, op, (oneKey) =>
            {
                Parallel.ForEach(IBDs[oneKey], (oneIBD) =>
                {
                    int sIndex = sitDict.Get_SiteIndex(oneIBD.Start);
                    int eIndex = sitDict.Get_SiteIndex(oneIBD.End);
                    for (int i = sIndex; i <= eIndex; i++)
                    {
                        Interlocked.Increment(ref cnts[i]);
                    }
                }
                );


            }
            );

            List<int> result = new List<int>();

            foreach (int oneCnt in cnts)
            {
                result.Add(oneCnt);
            }

            utl.listToFile(result, outPath);

            return result;
        }
    }
}
