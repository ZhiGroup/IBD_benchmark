using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace IBD_BM
{
    class Relatedness
    {



        /// <summary>
        /// This is a wraper for simple demo
        /// You may direct use Run method
        /// Need a supported data input format.
        /// You may add ore formats in the Loader module.
        /// 
        /// This module will output 4 numbers in console
        /// As # of individuals in:
        /// 4th degree, 3rd degree, 2nd degree, and 1st degree.
        /// 
        /// </summary>
        public static void singleRun()
        {
            string gMap_Path = "genetic map file";

            utl.GenMapV3 gMap = new utl.GenMapV3(gMap_Path, 1, 3);

            string inPath = "IBD output file";

            //vvvv change the data format here vvvvv
            Loader.dataType inType = Loader.dataType.GroudTruth;

            Run(inType, inPath, gMap);


        }

        /// <summary>
        /// actual core method
        /// </summary>
        /// <param name="inType"></param>
        /// <param name="path"></param>
        /// <param name="gMap"></param>
        /// <param name="IBD_Thd"></param>
        /// <returns></returns>
        static List<int> Run(Loader.dataType inType, string path, utl.GenMapV3 gMap, double IBD_Thd = 2)
        {
            Dictionary<long, double> Dict_indvPair_Len = new Dictionary<long, double>();

            string line;
            string[] parts;
            Loader.LoaderConfig cfg = new Loader.LoaderConfig(inType);
            //utl.filePercentageProgress prog = new utl.filePercentageProgress(path,10);
            StreamReader sr = new StreamReader(path);
            long key;
            double len;
            Loader.IBD_Phy_Start_End temIBD;

            if (inType == Loader.dataType.GroudTruth || inType == Loader.dataType.TPBWT)
            {
                line = sr.ReadLine();
            }




            while ((line = sr.ReadLine()) != null)
            {
                //prog.AddLine(line);
                parts = line.Split(cfg.Delimiter);
                key = Loader.MakeLongKey_noHapID(inType, parts);
                temIBD = new Loader.IBD_Phy_Start_End(parts, cfg, inType);
                len = Loader.Read_gLen(inType, cfg, parts);

                if (len < IBD_Thd)
                { continue; }

                if (Dict_indvPair_Len.ContainsKey(key) == false)
                {
                    Dict_indvPair_Len.Add(key, 0);
                }

                Dict_indvPair_Len[key] += len;


            }

            sr.Close();
            List<double> bins = new List<double>();
            bins.Add(gMap.maxGen() * 4 * (1 / Math.Pow(2, 11.0 / 2.0)));
            bins.Add(gMap.maxGen() * 4 * (1 / Math.Pow(2, 9.0 / 2.0)));
            bins.Add(gMap.maxGen() * 4 * (1 / Math.Pow(2, 7.0 / 2.0)));
            bins.Add(gMap.maxGen() * 4 * (1 / Math.Pow(2, 5.0 / 2.0)));
            Console.WriteLine(gMap.maxGen());
            Console.WriteLine(path);
            foreach (double one in bins)
            {
                Console.Write(one + "\t");
            }
            Console.WriteLine();

            List<int> cnts = utl.Distribution_GivenBins(Dict_indvPair_Len.Values.ToList(), bins);


            foreach (int oneCnt in cnts)
            {
                Console.Write(oneCnt + "\t");
            }
            Console.WriteLine();

            return cnts;
        }
    }
}
