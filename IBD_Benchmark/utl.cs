/*
Author: Kecong Tang(Benny)
Supporting Utility module, contains commonly used and reusealbe functions.
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Threading;

namespace IBD_BM
{
    class utl
    {
        public class GenMapV3
        {
            int indexCol = 0;
            int rateCol = 0;

            /// <summary>
            /// 
            /// </summary>
            /// <param name="mapPath"></param>
            /// <param name="phyCol_Index">0 based,the column has physical location</param>
            /// <param name="genCol_Index">0 based,the column has genetic rates</param>
            /// <param name="delimiter"></param>
            public GenMapV3(string mapPath, int phyCol_Index, int genCol_Index, char delimiter = '\t')
            {
                indexCol = phyCol_Index;
                rateCol = genCol_Index;
                loadGenMap(mapPath, delimiter);
            }



            Dictionary<int, double> rateDic = new Dictionary<int, double>();

            List<int> sorted_Keys = new List<int>();

            void loadGenMap(string dicFileName, char delimiter)
            {

                Console.WriteLine("Loading G-Map: " + dicFileName + "...");
                StreamReader sr = new StreamReader(dicFileName);
                string line;
                string[] parts;
                int index;
                double rate = 0, prevRate = double.MinValue;
                while ((line = sr.ReadLine()) != null
                    && line.Contains("Position") == true)
                { }

                do
                {
                    parts = line.Split(delimiter);

                    index = Convert.ToInt32(parts[indexCol]);
                    try
                    {
                        rate = Double.Parse(parts[rateCol], System.Globalization.NumberStyles.Any);
                    }
                    catch
                    {
                        rate = 0;
                        Console.WriteLine("Bug Decimal Convertion: " + line);
                    }



                    if (rate <= prevRate)
                    {
                        continue;
                    }
                    prevRate = rate;

                    rateDic.Add(index, rate);


                } while ((line = sr.ReadLine()) != null);
                sr.Close();


                List<int> keys = rateDic.Keys.ToList();
                keys.Sort();

                rate = 0;
                prevRate = double.MinValue;
                Dictionary<int, double> rateDicTem = new Dictionary<int, double>();
                foreach (int one in keys)
                {
                    if (rateDic[one] <= prevRate)
                    {
                        continue;
                    }
                    prevRate = rateDic[one];
                    rateDicTem.Add(one, rateDic[one]);
                }

                rateDic.Clear();

                rateDic = new Dictionary<int, double>(rateDicTem);

                rateDicTem.Clear();

                sorted_Keys = rateDic.Keys.ToList();
                sorted_Keys.Sort();

            }

            public int minPhy()
            {
                return rateDic.Keys.Min();
            }
            public int maxPhy()
            {
                return rateDic.Keys.Max();
            }

            public int LengthPhy()
            {
                return maxPhy() - minPhy();
            }

            public double LengthGen()
            {
                return maxGen() - minGen();
            }

            public double minGen()
            {
                return rateDic[minPhy()];
            }
            public double maxGen()
            {
                return rateDic[maxPhy()];
            }

            public double getGenLoc(int x)
            {
                if (rateDic.ContainsKey(x))
                {
                    return rateDic[x];
                }

                int prev = 0, next = 0;
                //int prevA = 0, nextA = 0;
                if (x < sorted_Keys.First())
                {
                    return 0;
                }
                else if (x > sorted_Keys.Last())
                {

                    return rateDic[sorted_Keys.Last()];
                }
                else
                {
                    int index = sorted_Keys.BinarySearch(x);

                    prev = sorted_Keys[~index - 1];
                    next = sorted_Keys[~index];
                }

                double result = rateDic[prev] + (x - prev) * (rateDic[next] - rateDic[prev]) / (next - prev);

                //rateDic.TryAdd(x, result);

                return result;


            }

            public double getGenLoc(string Pos)
            {
                int x = Convert.ToInt32(Pos);
                return getGenLoc(x);
            }

            public double getGenDistrance(int startPos, int endPos)
            {
                return getGenLoc(endPos) - getGenLoc(startPos);
            }

            public double getGenDistrance(string startPos, string endPos)
            {
                return getGenLoc(Convert.ToInt32(endPos)) - getGenLoc(Convert.ToInt32(startPos));
            }

        }

        public class siteDict
        {
            List<int> phyKeys = new List<int>();
            List<int> siteIndex_To_Phy = new List<int>();
            Dictionary<int, int> phy_To_SiteIndex = new Dictionary<int, int>();
            public siteDict(string vcfPath)
            {
                List<string> pos = get_All_POS(vcfPath);
                siteIndex_To_Phy = pos.Select(int.Parse).ToList();
                int cnt = 0;
                foreach (int one in siteIndex_To_Phy)
                {
                    phy_To_SiteIndex.Add(one, cnt);
                    cnt++;
                }
                phyKeys = new List<int>(phy_To_SiteIndex.Keys);
            }

            public int Count()
            {
                return phyKeys.Count();
            }
            public int Get_SiteIndex(int phy, bool startPhy = false)
            {
                if (phy_To_SiteIndex.ContainsKey(phy))
                {

                    return phy_To_SiteIndex[phy];
                }
                else
                {
                    int index = phyKeys.BinarySearch(phy);
                    if (startPhy)
                    {
                        if (~index - 1 < 0)
                        {
                            return phy_To_SiteIndex[phyKeys[0]];
                        }
                        else
                        {
                            return phy_To_SiteIndex[phyKeys[~index - 1]];
                        }
                    }
                    else
                    {
                        if (~index + 1 >= phyKeys.Count())
                        {
                            return phy_To_SiteIndex[phyKeys[phyKeys.Count() - 1]];
                        }
                        else
                        {
                            return phy_To_SiteIndex[phyKeys[~index + 1]];
                        }

                    }

                }
            }

            public int Get_Phy(int siteIndex)
            {
                return siteIndex_To_Phy[siteIndex];
            }

            static List<string> get_All_POS(string vcfPath)
            {
                List<string> pos = new List<string>();
                string line;
                string[] parts;
                StreamReader sr = new StreamReader(vcfPath);

                while ((line = sr.ReadLine()) != null && line.StartsWith("#"))
                {
                    continue;
                }

                do
                {
                  
                    parts = line.Substring(0, 30).Split('\t');
                    pos.Add(parts[1]);
                } while ((line = sr.ReadLine()) != null);
                sr.Close();
                return pos;
            }
        }

        /// <summary>
        /// check how much Range 1 is covered by Range 2
        /// </summary>
        /// <param name="r1_left"></param>
        /// <param name="r1_right"></param>
        /// <param name="r2_left"></param>
        /// <param name="r2_right"></param>
        /// <returns></returns>
        public static double rangeCovered_Percentage(double r1_left, double r1_right,
            double r2_left, double r2_right)
        {
            double len = Math.Min(r1_right, r2_right) - Math.Max(r1_left, r2_left);
            if (len < 0)
            {
                return 0;
            }

            return len / (r1_right - r1_left);
        }

        public static void listToFile<T>(List<T> data, string filePath)
        {
            Console.WriteLine("Writing " + filePath);
            StreamWriter sw = new StreamWriter(filePath);
            foreach (var one in data)
            {
                sw.WriteLine(one);
            }
            sw.Close();
        }


        /// <summary>
        /// given a list of bin
        /// bins must be increment sorted
        /// </summary>
        /// <param name="nums"></param>
        /// <param name="bins"></param>
        /// <returns></returns>
        public static List<int> Distribution_GivenBins(List<double> nums, List<double> bins)
        {
            //List<int> res = new List<int>(bins.Count());
            int[] res = new int[bins.Count()];
            Parallel.ForEach(nums, (oneNum) =>
            //foreach (double oneNum in nums)
            {
                for (int BI = bins.Count() - 1; BI >= 0; BI--)
                {
                    if (oneNum > bins[BI])
                    {
                        Interlocked.Increment(ref res[BI]);
                        //res[BI]++;
                        break;
                    }
                }


            });

            return res.ToList();
        }
    }
}
