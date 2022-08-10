using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace DownSample
{
    class MapGroundTruth
    {
        class L_R
        {
            public int L_Val = 0;
            public int R_Val = 0;
            public L_R(int left, int right)
            {
                L_Val = left;
                R_Val = right;
            }

        }

        /// <summary>
        /// give a smaller set of markers(arrVCF) and a bigger set of markers(seqVCF)
        /// make a ground truth file for the smaller set of marker.
        /// bigger set must fully include the smaller set.
        /// </summary>
        /// <param name="arrVCF_Path">low density vcf, no necessary to be array data</param>
        /// <param name="seqVCF_Path">high density vcf, no necessary to be sequencing</param>
        /// <param name="gtPath"></param>
        public static void run(string arrVCF_Path, string seqVCF_Path, string gtPath, string outPath)
        {


            Dictionary<int, L_R> Seq_To_Arr = new Dictionary<int, L_R>();
            HashSet<int> arrHash = new HashSet<int>();

            List<int> arrPos = new List<int>();
            List<int> seqPos = new List<int>();

            arrPos = utl.get_All_POS(arrVCF_Path).Select(int.Parse).ToList();
            seqPos = utl.get_All_POS(seqVCF_Path).Select(int.Parse).ToList();
            arrHash = new HashSet<int>(arrPos);
            arrPos.Sort();

            int index;

            int prev;
            int next;
            foreach (int oneSeq in seqPos)
            {
                if (arrHash.Contains(oneSeq))
                {
                    continue;
                }

                if (oneSeq < arrPos.First())
                {
                    Seq_To_Arr.Add(oneSeq, new L_R(arrPos.First(), arrPos.First()));
                    continue;
                }

                if (oneSeq > arrPos.Last())
                {
                    Seq_To_Arr.Add(oneSeq, new L_R(arrPos.Last(), arrPos.Last()));
                    continue;
                }


                index = arrPos.BinarySearch(oneSeq);

                prev = arrPos[~index - 1];
                next = arrPos[~index];

                Seq_To_Arr.Add(oneSeq, new L_R(prev, next));

            }

            //#individual_1_id,individual_1_haplotype_id,individual_2_id,individual_2_haplotype_id,chromosome_id,true_ibd_physical_position_start,true_ibd_physical_position_end,genetic_length
            //342,1,7,0,20,66785,353674,1.034238
            //650,0,14,1,20,66785,353674,1.034238


            string line;
            string[] parts;
            StreamReader sr = new StreamReader(gtPath);
            StreamWriter sw = new StreamWriter(outPath);
            sw.NewLine = "\n";
            line = sr.ReadLine();
            sw.WriteLine(line);
            int sPos, ePos;
            while ((line = sr.ReadLine()) != null)
            {
                parts = line.Split(',');
                sPos = Convert.ToInt32(parts[5]);
                ePos = Convert.ToInt32(parts[6]);
                if (arrHash.Contains(sPos) == false)
                {
                    sPos = Seq_To_Arr[sPos].R_Val;
                }

                if (arrHash.Contains(ePos) == false)
                {
                    ePos = Seq_To_Arr[ePos].L_Val;
                }
                for (int i = 0; i < 5; i++)
                {
                    sw.Write(parts[i] + ",");
                }
                sw.WriteLine(sPos + "," + ePos + "," + parts[7]);
            }



            sw.Close();
            sr.Close();


        }


    }
}
