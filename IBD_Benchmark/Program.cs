/*
Author: Kecong Tang(Benny)
Program entrance, for clear demonstration, only a single example is provied.

Three major entrances availables in this package:
1. Call GroupCaller.run() to run power and accuracy measures.
2. Call Relatedness.singleRun() for relatedness analysis.
3. Call siteCoverage.singleRun() for IBD distribution.

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

        public static bool parallel_Compute = true;

        public static bool phasingErrorData = false;
        public static double minBin = 2;
        public static double maxBin = 7;
        public static double binLen = 1;

        public static string tool_Name = "";
        public static Loader.dataType tool_Type;

        public static string reported_Path = "";
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
            Console.WriteLine("IBD Benchmark Tool v1.01.");

            loadConfig(args);

            GroupCaller.run();

        }

        public static void loadConfig(string[] args)
        {
            if (args.Count() < 1 )
            {
                if (File.Exists(configFilePath) == false)
                {
                    Console.WriteLine("Missing Configuration File Path.");
                    Environment.Exit(1);
                }

            }
            else
            {
                if (string.IsNullOrWhiteSpace(args[0]) == false)
                {
                    configFilePath = args[0];
                }
                else
                {
                    Console.WriteLine("Missing Configuration File Path.");
                    Environment.Exit(1);
                }
            }
     

            Console.WriteLine("Loading "+configFilePath+" ...");
            string[] lines = File.ReadAllLines(configFilePath);
            string[] parts;
            foreach (string line in lines)
            {
                parts = line.Split('\t');


                if (line.StartsWith("phasingErrorData"))
                {
                    phasingErrorData = boolCheck_Read("phasingErrorData", parts);
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


                if (line.StartsWith("reported_Path"))
                {

                    reported_Path = pathCheck_Read("reported_Path", parts);
                    continue;
                }

                if (line.StartsWith("tool_Name"))
                {
                    parts = line.Split('\t');
                    if (parts.Count() < 2)
                    {
                        Console.WriteLine("tool_Name not given, please provide tool name and and make sure using Tab delimiter.");
                        Environment.Exit(1);
                    }

                    tool_Name = parts[1].Trim();

                    switch (tool_Name.ToLower())
                    {
                        case "fastsmc":
                            tool_Type = Loader.dataType.FastSMC;
                            break;
                        case "hap-ibd":
                            tool_Type = Loader.dataType.HapIBD;
                            break;
                        case "ilash":
                            tool_Type = Loader.dataType.iLash;
                            break;
                        case "rapid":
                            tool_Type = Loader.dataType.RaPID;
                            break;
                        case "tpbwt":
                            tool_Type = Loader.dataType.TPBWT;
                            break;
                        default:
                            Console.WriteLine("No Supported Tool Type["+tool_Name+"]!");
                            Console.WriteLine("Currently supported formats: FastSMC,hap-IBD,iLash,RaPID and TPBWT, add/modify the paser in Loader.cs for more formats");
                            Environment.Exit(1);
                            break;

                    }
                    Console.WriteLine("IBD Tool:\t" + tool_Name);
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
                    out_Path = dirCheck_Read("out_Path", parts);
                    continue;
                }





            }

            if (String.IsNullOrWhiteSpace(out_Path))
            {
                out_Path = reported_Path + ".IBD_BM_Result.txt";
            }
            Console.WriteLine(configFilePath + " Loaded");



        }


        public static bool boolCheck_Read(string name, string[] parts)
        {


            if (parts.Count() != 2 || String.IsNullOrWhiteSpace(parts[1]))
            {
                Console.WriteLine("Path Error [" + name + "]. Check path and make sure using Tab delimiter.");
                Console.ReadKey();
                Environment.Exit(1);
            }

            Console.WriteLine(name + ":\t" + parts[1]);

            if (parts[1].ToLower() == "true")
            {
                return true;
            }
            else
            {
                return false;
            }
    
        }

        public static string dirCheck_Read(string name, string[] parts)
        {


            if (parts.Count() != 2 || Directory.Exists(Path.GetDirectoryName( parts[1])) == false || String.IsNullOrWhiteSpace(parts[1]))
            {
                Console.WriteLine("Path Error [" + name + "]. Check path and make sure using Tab delimiter.");
                Console.ReadKey();
                Environment.Exit(1);
            }

            Console.WriteLine(name + ":\t" + parts[1]);
            return parts[1];
        }


        public static string pathCheck_Read(string name,string[] parts)
        {
            

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
