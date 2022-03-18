/*
Author: Kecong Tang(Benny)
Program entrance, for clear demonstration, only a single example is provied.
The defalut setting does parallel hard drive reading, may be disable by modifying Config.txt.
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace IBD_BM
{

    class Program
    {
        #region Configurations could be modified in Config.txt
        public static bool parallel_Load = true;
        public static bool parallel_Compute = true;

        public static bool phasingErrorData = false;
        public static double minBin = 2;
        public static double maxBin = 7;
        public static double binLen = 1;

        public static string RP_Path = "";
        public static string HI_Path = "";
        public static string TP_Path = "";
        public static string IL_Path = "";
        public static string FS_Path = "";
        public static string gMap_Path = "";
        public static string vcf_Path = "";
        public static string GT_Path = "";
        public static string out_Path = "";
        public static int gMap_PositionCol_Index_ZeroBased = 1;
        public static int gMap_MapCol_Index_ZeroBased = 3;
        #endregion

        public static string configFilePath = "Config.txt";
        
        static void Main(string[] args)
        {
            loadConfig();

            GroupCaller.run();

        }

        public static void loadConfig()
        {
      
            Console.WriteLine("Loading "+configFilePath+" ...");
            string[] lines = File.ReadAllLines(configFilePath);
            string[] parts;
            foreach (string line in lines)
            {
                parts = line.Split('\t');


                if (line.StartsWith("parallel_Load"))
                {
                    if (parts[1].ToLower() == "true")
                    {
                        parallel_Load = true;
                    }
                    else if (parts[1].ToLower() == "false")
                    {
                        parallel_Load = false;
                    }
                    Console.WriteLine(line);
                    continue;
                }

                if (line.StartsWith("phasingErrorData"))
                {
                    if (parts[1].ToLower() == "true")
                    {
                        phasingErrorData = true;
                    }
                    else if (parts[1].ToLower() == "false")
                    {
                        phasingErrorData = false;
                    }
                    Console.WriteLine(line);
                    continue;
                }


                if (line.StartsWith("gMap_PositionCol_Index_ZeroBased"))
                {
                    gMap_PositionCol_Index_ZeroBased = Convert.ToInt32(parts[1]);
                    Console.WriteLine(line);
                    continue;
                }

                if (line.StartsWith("gMap_MapCol_Index_ZeroBased"))
                {
                    gMap_MapCol_Index_ZeroBased = Convert.ToInt32(parts[1]);
                    Console.WriteLine(line);
                    continue;
                }

                if (line.StartsWith("minBin"))
                {
                    minBin = Convert.ToInt32(parts[1]);
                    Console.WriteLine(line);
                    continue;
                }

                if (line.StartsWith("maxBin"))
                {
                    maxBin = Convert.ToInt32(parts[1]);
                    Console.WriteLine(line);
                    continue;
                }

                if (line.StartsWith("binLen"))
                {
                    binLen = Convert.ToInt32(parts[1]);
                    Console.WriteLine(line);
                    continue;
                }


                if (line.StartsWith("FS_Path"))
                {

                    FS_Path = pathCheck_Read("FS_Path", parts);
                    continue;
                }
                if (line.StartsWith("HI_Path"))
                {

                    HI_Path = pathCheck_Read("HI_Path", parts);
                    continue;
                }
                if (line.StartsWith("IL_Path"))
                {

                    IL_Path = pathCheck_Read("IL_Path", parts);
                    continue;
                }
                if (line.StartsWith("RP_Path"))
                {

                    RP_Path = pathCheck_Read("RP_Path", parts);
                    continue;
                }
                if (line.StartsWith("TP_Path"))
                {

                    TP_Path = pathCheck_Read("TP_Path", parts);
                    continue;
                }

                if (line.StartsWith("gMap_Path"))
                {
                    gMap_Path = pathCheck_Read("gMap_Path", parts);
                    continue;
                }
                if (line.StartsWith("vcf_Path"))
                {
                    vcf_Path = pathCheck_Read("vcf_Path", parts);
                    continue;
                }
                if (line.StartsWith("GT_Path"))
                {
                    GT_Path = pathCheck_Read("GT_Path", parts);
                    continue;
                }
                if (line.StartsWith("out_Path"))
                {
                    out_Path = pathCheck_Read("out_Path", parts);
                    continue;
                }

            }


            Console.WriteLine(configFilePath + " Loaded");



        }

        public static string pathCheck_Read(string name,string[] parts)
        {
            if (name == "out_Path")
            {
                if (parts.Count() != 2 || String.IsNullOrWhiteSpace(parts[1]))
                {
                    Console.WriteLine("Path Error [" + name + "]. Check path and make sure using Tab delimiter.");
                    Console.ReadKey();
                    Environment.Exit(1);
                }
                Console.WriteLine(name + ":\t" + parts[1]);
                return parts[1];
            }
            

            if (parts.Count()!=2|| File.Exists(parts[1]) == false || String.IsNullOrWhiteSpace(parts[1]))
            {
                Console.WriteLine("Path Error [" + name + "]. Check path and make sure using Tab delimiter.");
                Console.ReadKey();
                Environment.Exit(1);
            }

            Console.WriteLine(name + ":\t" + parts[1]);
            return parts[1];
        }

    }
}
