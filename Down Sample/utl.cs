using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace DownSample
{
    class utl
    {
        public static int getPOS_Int(string line)
        {
            string pos = getPOS(line);
            return Convert.ToInt32(pos);

        }

        public static string getPOS(string line)
        {
            //#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NWD112649	NWD278543	
            //         1   2   3   4   5     6         7     8         9

            int posStart = 0;
            int posEnd = 0;
            int tabCnt = 0;
            for (int i = 0; i < line.Length; i++)
            {
                if (line[i] == '\t')
                {
                    tabCnt++;
                    if (tabCnt == 1)
                    {
                        posStart = i + 1;
                    }
                    if (tabCnt == 2)
                    {
                        posEnd = i;
                        break;
                    }
                }



            }
            return line.Substring(posStart, posEnd - posStart);
        }
        public static List<string> get_All_POS(string vcfPath)
        {
            List<string> pos = new List<string>();
            string line;
            StreamReader sr = new StreamReader(vcfPath);
            while ((line = sr.ReadLine()) != null && line.StartsWith("#"))
            {
                continue;
            }

            do
            {
                pos.Add(getPOS(line));
            } while ((line = sr.ReadLine()) != null);
            sr.Close();
            return pos;
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



        public static int lineCount(string path, string skipHeader = "")
        {
            Console.WriteLine("Counting Lines...\t" + path);
            StreamReader sr = new StreamReader(path);
            string line;
            int cnt = 0;
            if (String.IsNullOrEmpty(skipHeader))
            {
                while ((line = sr.ReadLine()) != null)
                {
                    cnt++;
                }
            }
            else
            {
                while ((line = sr.ReadLine()) != null && line.StartsWith(skipHeader))
                {
                    continue;
                }
                cnt++;
                while ((line = sr.ReadLine()) != null)
                {
                    cnt++;
                }

            }
            sr.Close();
            Console.WriteLine(cnt + "\tlines");
            return cnt;
        }


        public class GenMapV3
        {
            int indexCol = 0;
            int rateCol = 0;
            Dictionary<int, double> rateDic = new Dictionary<int, double>();

            List<int> sorted_Keys = new List<int>();
            public List<int> GetPosInt()
            {
                return sorted_Keys;
            }
            public List<string> GetPosStr()
            {
                return (from i in sorted_Keys select i.ToString()).ToList();
            }

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
        /// <summary>
        /// the first index is column index, the second is the row index
        /// </summary>
        /// <param name="path"></param>
        /// <param name="nthCols">all columns that need to collect 0 based!!</param>
        /// <param name="delimiter"></param>
        /// <returns></returns>
        public static List<string> fileToList(string path, int rowsToSkip = 0, bool ignoreNullString = true)
        {
            List<string> results = new List<string>();

            StreamReader sr = new StreamReader(path);
            string line;

            for (int i = 0; i < rowsToSkip; i++)
            {
                line = sr.ReadLine();
            }
            while ((line = sr.ReadLine()) != null)
            {
                if (ignoreNullString == true
                    && String.IsNullOrWhiteSpace(line))
                {
                    continue;
                }
                results.Add(line);

            }

            return results;
        }

        public static int MinorAlleleCnt(string[] parts)
        {
            int oneCnt = 0;
            int zeroCnt = 0;
            for (int i = 9; i < parts.Count(); i++)
            {
                if (parts[i][0] == '1')
                {
                    oneCnt++;
                }
                else if (parts[i][0] == '0')
                {
                    zeroCnt++;
                }

                if (parts[i][2] == '1')
                {
                    oneCnt++;
                }
                else if (parts[i][2] == '0')
                {
                    zeroCnt++;
                }

            }

            return Math.Min(oneCnt, zeroCnt);

        }


        /// <summary>
        /// divide a list into n chunks.
        /// there many be a little more item in the last chunk
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="input"></param>
        /// <param name="nChunks"></param>
        /// <returns></returns>
        public static List<List<T>> DivideList<T>(List<T> input, int nChunks)
        {

            List<List<T>> res = new List<List<T>>();

            if (nChunks == 0)
            { return res; }

            if (nChunks > input.Count())
            {
                Console.WriteLine("DivideList Wrong Input! input.Count() > nChunks");

                foreach (T one in input)
                {
                    List<T> tem = new List<T>();
                    tem.Add(one);
                    res.Add(tem);
                }

                return res;

            }

            int wSize = input.Count() / nChunks;
            for (int i = 0; i < nChunks; i++)
            {
                res.Add(input.GetRange(i * wSize, wSize));
            }

            //tail
            for (int i = wSize * nChunks; i < input.Count(); i++)
            {
                res.Last().Add(input[i]);
            }

            return res;
        }
    }
}
