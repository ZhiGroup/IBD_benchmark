using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Collections;

namespace makeGT_Err
{
    class Program
    {
        static void Main(string[] args)
        {

            string inPath = "E:\\tem\\dummy.vcf";

            List<double> rates = new List<double>();
            rates.Add(0.003);
            rates.Add(0.002);
            rates.Add(0.001);
            //add rate as needed

            makeVCF(inPath, rates);

        }


        /// <summary>
        /// create pannel with error rates
        /// </summary>
        /// <param name="inPath"></param>
        public static void makeVCF(string inPath, List<double> rates)
        {

            int nIndv = utl.get_nIndv(inPath);
            int nPos = utl.get_All_POS(inPath).Count();


            List<int> nToTake = new List<int>();
            foreach (double one in rates)
            {
                nToTake.Add(Convert.ToInt32(one / 2 * nPos));
            }

            Console.WriteLine("Making Error Tables...");
            List<List<int>> errTable = makeErrorIndex(rates.Max()/2, nIndv * 2, nPos);


            Console.WriteLine("Making Index Tables...");
            List<List<HashSet<int>>> errIndexTables = new List<List<HashSet<int>>>();
            ///make hash tables
            for (int r = 0; r < rates.Count(); r++)
            {
                errIndexTables.Add(new List<HashSet<int>>());

                for (int i = 0; i < nIndv * 2; i++)
                {
                    errIndexTables[r].Add(new HashSet<int>(errTable[i].Take(nToTake[r])));
                }
            }

            List<StreamWriter> SWs = new List<StreamWriter>();
            for (int i = 0; i < errIndexTables.Count(); i++)
            {
                string name = inPath.Replace(".vcf", ".e") + rates[i].ToString() + ".vcf";
                SWs.Add(new StreamWriter(name));
                SWs[i].NewLine = "\n";
            }


            StreamReader sr = new StreamReader(inPath);
            string line;
            string[] parts;
            while ((line = sr.ReadLine()) != null && line.StartsWith("##"))
            {
                for (int i = 0; i < errIndexTables.Count(); i++)
                {
                    SWs[i].WriteLine(line);
                }
            }
            for (int i = 0; i < errIndexTables.Count(); i++)
            {
                SWs[i].WriteLine(line);
            }



            int rowCnt = 0;
            while ((line = sr.ReadLine()) != null)
            {
                if (rowCnt % 10000 == 0)
                {
                    Console.WriteLine(DateTime.Now + " " + rowCnt + "/" + nPos);
                }
                parts = line.Split('\t');
                Parallel.For(0, errIndexTables.Count(), (r) =>
                //for (int r = 0; r < errIndexTables.Count(); r++)
                {
                    StringBuilder sb = new StringBuilder();
                    int hapID = 0;
                    for (int c = 0; c < parts.Count(); c++)
                    {
                        if (c < 9)
                        {
                            sb.Append(parts[c] + "\t");
                        }
                        else
                        {
                            //indv ID to hap ID
                            hapID = (c - 9) * 2;
                            if (errIndexTables[r][hapID].Contains(rowCnt))
                            {
                                sb.Append(ZeroOneShift(parts[c][0]));
                            }
                            else
                            {
                                sb.Append(parts[c][0]);
                            }
                            sb.Append("|");

                            hapID = hapID + 1;
                            if (errIndexTables[r][hapID].Contains(rowCnt))
                            {
                                sb.Append(ZeroOneShift(parts[c][2]));
                            }
                            else
                            {
                                sb.Append(parts[c][2]);
                            }
                            sb.Append("\t");
                        }
                    }
                    sb.Length--;
                    SWs[r].WriteLine(sb.ToString());
                });
                rowCnt++;
            }

            sr.Close();
            for (int i = 0; i < errIndexTables.Count(); i++)
            {
                SWs[i].Close();
            }
        }

        public static char ZeroOneShift(char c)
        {
            if (c == '1')
            {
                return '0';
            }
            else
            {
                return '1';
            }


        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="rate_PerRow"></param>
        /// <param name="nRow">indv x 2</param>
        /// <param name="nCell_PerRow">pos</param>
        /// <returns>inner List<int> is a hap </returns>
        public static List<List<int>> makeErrorIndex(double rate_PerRow, int nHap, int nSite)
        {
            List<List<int>> result = new List<List<int>>();
            int nSelect = Convert.ToInt32(rate_PerRow * nSite);
            for (int i = 0; i < nHap; i++)
            {
                result.Add(new List<int>());
            }

            int cnt = 0;
            ParallelOptions op = new ParallelOptions();
            op.MaxDegreeOfParallelism = 10;
            Parallel.For(0, nHap, op, (i) =>
            {
                int tID = System.Threading.Thread.CurrentThread.ManagedThreadId;
                result[i].AddRange(randomSelect_Fast(nSelect, nSite, tID + i));
                cnt++;
                if (cnt % 1000 == 0)
                {
                    Console.WriteLine(DateTime.Now + " " + cnt + "/" + nHap);
                }
            });

            return result;
        }

        /// <summary>
        /// random select nSelect from nTotal
        /// </summary>
        /// <param name="nSelect"></param>
        /// <param name="nTotal"></param>
        /// <returns></returns>
        public static List<int> randomSelect(int nSelect, int nTotal, int extraSeed = 1)
        {
            List<int> result = new List<int>();
            Random rnd = new Random((int)(DateTime.Now.Ticks + extraSeed));
            List<double> nums = new List<double>();
            Dictionary<double, int> val_To_Index = new Dictionary<double, int>();
            double oneNum = 0;

            int cnt = 0;
            while (nums.Count() < nTotal)
            {
                oneNum = rnd.NextDouble();
                if (val_To_Index.ContainsKey(oneNum))
                { continue; }

                val_To_Index.Add(oneNum, cnt);
                nums.Add(oneNum);
                cnt++;
            }


            nums.Sort();

            for (int i = 0; i < nSelect; i++)
            {
                result.Add(val_To_Index[nums[i]]);
            }

            return result;
        }


        /// <summary>
        /// random select nSelect from nTotal
        /// </summary>
        /// <param name="nSelect"></param>
        /// <param name="nTotal"></param>
        /// <returns></returns>
        public static List<int> randomSelect_Fast(int nSelect, int nTotal, int extraSeed = 1)
        {
            List<int> result = new List<int>();
            Random rnd = new Random((int)(DateTime.Now.Ticks + extraSeed));
            HashSet<int> chosen = new HashSet<int>();

            int oneNum;
            while (chosen.Count() < nSelect)
            {
                oneNum = rnd.Next(nTotal);
                chosen.Add(oneNum);
            }
            return chosen.ToList();
        }

    }
}
