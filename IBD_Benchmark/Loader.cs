/*
Author: Kecong Tang(Benny)
Loading unit, contains parsing methods for ground truth files and result files from IBD tools.
Parser could to be modified for different formats.
*/
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace IBD_Benchmark
{
    class Loader
    {
        public static int ID_Key_Pad_Num = 1000000;
        /// <summary>
        /// ID padding: 976	3285	1	0//ID1 ID2 hap1 hap2
        /// (1000000+ID1)x10_hap1_00_(1000000+ID2)x10_hap2
        /// 10000000              00  10000000
        /// 100,000,000,010,000,000
        /// </summary>
        public enum dataType { GroudTruth, RaPID, Ger2, HapIBD, TPBWT, iLash, FastSMC };
        public class LoaderConfig
        {
            public int ID1_idx = 0;
            public int ID2_idx = 0;
            public int Hap1_idx = 0;
            public int Hap2_idx = 0;
            public int PhyStart_idx = 0;
            public int PhyEnd_idx = 0;
            public int Len_idx = -1;
            public bool SkipFirstRow = false;
            public char Delimiter = '\t';
            public LoaderConfig(dataType inType)
            {
                switch (inType)
                {
                    case dataType.GroudTruth:
                        //#individual_1_id,individual_1_haplotype_id,individual_2_id,individual_2_haplotype_id,chromosome_id,true_ibd_physical_position_start,true_ibd_physical_position_end,genetic_length
                        //3985,1,3819,0,20,62062477,62948300,1.608621429017873
                        ID1_idx = 0;
                        ID2_idx = 2;
                        Hap1_idx = 1;
                        Hap2_idx = 3;
                        PhyStart_idx = 5;
                        PhyEnd_idx = 6;
                        Len_idx = 7;
                        SkipFirstRow = true;
                        Delimiter = ',';
                        break;
                    case dataType.RaPID:
                        //21	976	3285	1	0	24412730	27809303	5.42485	32175	36539
                        ID1_idx = 1;
                        ID2_idx = 2;
                        Hap1_idx = 3;
                        Hap2_idx = 4;
                        PhyStart_idx = 5;
                        PhyEnd_idx = 6;
                        Len_idx = 7;
                        break;
                    case dataType.TPBWT:
                        //21 20407 15508510 5.0415359279141105 1181 0 2873 0 0 1733 0.0
                        //0    1      2                3         4  5  6   7
                        //^^^^^old format^^^^^
                        //vvvvvnew formatvvvvv
                        //0,1,         2    3      4     5   6              7   8            9      10     11
                        //,chromosome,end,end_bp,end_cm,id1,id1_haplotype,id2,id2_haplotype,start,start_bp,start_cm
                        //0,20,785,1858299,5.07728238109226,835,0,2746,1,0,565687,2.8046106882790003e-07

                        ID1_idx = 5;
                        ID2_idx = 7;
                        Hap1_idx = 6;
                        Hap2_idx = 8;
                        PhyStart_idx = 10;
                        PhyEnd_idx = 3;
                        Delimiter = ',';
                        SkipFirstRow = true;
                        break;
                    case dataType.HapIBD:
                        //1402	2	2583	2	21	774646	15379203	10.154
                        ID1_idx = 0;
                        ID2_idx = 2;
                        Hap1_idx = 1;
                        Hap2_idx = 3;
                        PhyStart_idx = 5;
                        PhyEnd_idx = 6;
                        Len_idx = 7;
                        break;
                    case dataType.Ger2:
                        //0 2540	0 3666	21	39322540 44203275	. .	6361	12.68	cM	0	1	0
                        ID1_idx = 1;
                        ID2_idx = 3;
                        Hap1_idx = 0;
                        Hap2_idx = 2;
                        PhyStart_idx = 4;
                        PhyEnd_idx = 5;
                        Len_idx = 10;
                        break;

                    case dataType.iLash:
                        //0	3599_1	0	3961_1	20	50209309	62948300	.	.	15.7624	1
                        ID1_idx = 1;
                        ID2_idx = 3;
                        PhyStart_idx = 5;
                        PhyEnd_idx = 6;
                        Len_idx = 9;
                        break;
                    case dataType.FastSMC:
                        //1815	1815	1	3163	3163	2	20	621161	2235750	6.013828	0.9912891	67.70425	24.99993
                        //0. First individual's family identifier
                        //1. First individual identifier
                        //2. First individual haplotype identifier (1 or 2)
                        //3. Second individual's family identifier
                        //4. Second individual identifier
                        //5. Second individual haplotype identifier (1 or 2)
                        //6. Chromosome number
                        //7. Starting position of the IBD segment (inclusive)
                        //8. Ending position of the IBD segment (inclusive)
                        //9. (optional) Length in centimorgans of IBD segment
                        //10. IBD score
                        //11. (optional) Average mean posterior age estimate of the IBD segment
                        //12. (optional) Average MAP age estimate of the IBD segment
                        ID1_idx = 1;
                        ID2_idx = 3;
                        Hap1_idx = 2;
                        Hap2_idx = 5;
                        PhyStart_idx = 7;
                        PhyEnd_idx = 8;
                        Len_idx = 9;
                        break;



                    default:
                        Console.WriteLine("BUG!!Wrong Type!");
                        break;

                }

            }
        }

        public static long MakeLongKey(dataType inType, string[] parts)
        {

            //now all set have A0 ID

            /// ID padding: 976	3285	1	0//ID1 ID2 hap1 hap2
            /// (1000000+ID1)x10_hap1_00_(1000000+ID2)x10_hap2
            /// 100,000,00              00  10000000
            /// 100,000,000,010,000,000
            LoaderConfig cfg = new LoaderConfig(inType);

            if (parts[cfg.ID1_idx].StartsWith("A"))
            {
                parts[cfg.ID1_idx] = parts[cfg.ID1_idx].Substring(1);
            }

            if (parts[cfg.ID2_idx].StartsWith("A"))
            {
                parts[cfg.ID2_idx] = parts[cfg.ID2_idx].Substring(1);
            }

            int tem1 = 0, tem2 = 0;
            int ID1 = 0, ID2 = 0, hap1 = 0, hap2 = 0;
            switch (inType)
            {

                case dataType.HapIBD:
                    //737	2	2948	2	21	40438035	43459098	9.272
                    ID1 = Convert.ToInt32(parts[0]);
                    hap1 = Convert.ToInt32(parts[1]) - 1;
                    ID2 = Convert.ToInt32(parts[2]);
                    hap2 = Convert.ToInt32(parts[3]) - 1;
                    break;

                case dataType.iLash:
                    //0	2700_1	0	3998_0	20	52477650	56479537	.	.	5.01981	1

                    ID1 = Convert.ToInt32(parts[1].Split('_')[0]);
                    hap1 = Convert.ToInt32(parts[1].Split('_')[1]);
                    ID2 = Convert.ToInt32(parts[3].Split('_')[0]);
                    hap2 = Convert.ToInt32(parts[3].Split('_')[1]);
                    break;
                case dataType.FastSMC:
                    //1815	1815	1	3163	3163	2	20	621161	2235750	6.013828	0.9912891	67.70425	24.99993
                    //0. First individual's family identifier
                    //1. First individual identifier
                    //2. First individual haplotype identifier (1 or 2)
                    //3. Second individual's family identifier
                    //4. Second individual identifier
                    //5. Second individual haplotype identifier (1 or 2)
                    //6. Chromosome number
                    //7. Starting position of the IBD segment (inclusive)
                    //8. Ending position of the IBD segment (inclusive)
                    //9. (optional) Length in centimorgans of IBD segment
                    //10. IBD score
                    //11. (optional) Average mean posterior age estimate of the IBD segment
                    //12. (optional) Average MAP age estimate of the IBD segment

                    //i patched indv ID 0 to A0 since they don't allow 0 as indv ID.
                    //need to convert it back

                    ID1 = Convert.ToInt32(parts[1]);
                    ID2 = Convert.ToInt32(parts[4]);
                    hap1 = Convert.ToInt32(parts[2]) - 1;
                    hap2 = Convert.ToInt32(parts[5]) - 1;
                    break;
                default:
                    //format of GroundTruth TPBWT and RaPID
                    //case dataType.GroudTruth:
                    //#individual_1_id,individual_1_haplotype_id,individual_2_id,individual_2_haplotype_id,chromosome_id,true_ibd_physical_position_start,true_ibd_physical_position_end,genetic_length
                    //3985,1,3819,0,20,62062477,62948300,1.608621429017873
                    //rapid
                    //20	A0	126	0	1	59653905	62405693	8.21924	15723	17196
                    ID1 = Convert.ToInt32(parts[cfg.ID1_idx]);
                    ID2 = Convert.ToInt32(parts[cfg.ID2_idx]);
                    hap1 = Convert.ToInt32(parts[cfg.Hap1_idx]);
                    hap2 = Convert.ToInt32(parts[cfg.Hap2_idx]);
                    break;
            }

            int tem;
            if (ID2 < ID1)
            {
                tem = ID2;
                ID2 = ID1;
                ID1 = tem;
                tem = hap2;
                hap2 = hap1;
                hap1 = tem;
            }
            long high = (ID_Key_Pad_Num + ID1) * 10 + hap1;
            long low = (ID_Key_Pad_Num + ID2) * 10 + hap2;
            return high * ID_Key_Pad_Num * 10000 + low;

        }

        public class IBD_Phy_Start_End
        {
            public int Start;
            public int End;
            public IBD_Phy_Start_End()
            { }
            public IBD_Phy_Start_End(int startPhy, int endPhy)
            {
                Start = startPhy;
                End = endPhy;
            }
            public IBD_Phy_Start_End(string[] parts, LoaderConfig cfg, dataType inType)
            {
                //ger
                //0 52	0 189	21	1733 5016422	. .	6279	5.01	MB	0	1	1
                //    tab    tab  tab           tab    tab  tab     tab   ......
                if (inType == dataType.Ger2)
                {
                    Start = Convert.ToInt32(parts[3].Split(' ')[0]);
                    End = Convert.ToInt32(parts[3].Split(' ')[1]);
                }
                else
                {
                    Start = Convert.ToInt32(parts[cfg.PhyStart_idx]);
                    End = Convert.ToInt32(parts[cfg.PhyEnd_idx]);
                }
            }

        }

        public static Dictionary<long, List<IBD_Phy_Start_End>> Load_IBD(dataType inType, string path, double minLen = 0)
        {
            Dictionary<long, List<IBD_Phy_Start_End>> holder = new Dictionary<long, List<IBD_Phy_Start_End>>();
            LoaderConfig cfg = new LoaderConfig(inType);
            string line;
            string[] parts;

            StreamReader sr = new StreamReader(path);

            if (cfg.SkipFirstRow == true)
            {
                line = sr.ReadLine();
            }
            int errorCnt = 0;
            long key;
            double minLenREP = 0;
            while ((line = sr.ReadLine()) != null)
            {

                parts = line.Split(cfg.Delimiter);
                if (inType == dataType.TPBWT)
                {
                    minLenREP = Convert.ToDouble(parts[4]) - Convert.ToDouble(parts[11]);
                }
                else
                {
                    minLenREP = Convert.ToDouble(parts[cfg.Len_idx]);
                }

                if (minLen != 0 && minLenREP < minLen)
                { continue; }
                try
                {
                    key = MakeLongKey(inType, parts);
                    if (!holder.ContainsKey(key))
                    {
                        holder.Add(key, new List<IBD_Phy_Start_End>());
                    }

                    holder[key].Add(new IBD_Phy_Start_End(parts, cfg, inType));
                }
                catch
                {
                    errorCnt++;
                }
            }
            sr.Close();

            return holder;
        }
    }
}
