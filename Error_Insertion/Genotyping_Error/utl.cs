using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;


namespace makeGT_Err
{
    class utl
    {
        public static int get_nIndv(string path)
        {
            StreamReader sr = new StreamReader(path);
            string line;
            string[] parts;

            List<string> indvTags = new List<string>();

            while ((line = sr.ReadLine()) != null && line.StartsWith("#CHROM") == false)
            { continue; }

            int cutPoint = indexOf_nTH_Char(line, 9, '\t');
            parts = line.Substring(cutPoint).Trim().Split('\t');

            sr.Close();
            return parts.Count();

        }

        public static List<string> get_All_POS(string vcfPath)
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

        public static int indexOf_nTH_Char(string orgString, int n, char c)
        {
            int result = -1;
            int cnt = 0;
            for (int i = 0; i < orgString.Length; i++)
            {
                if (orgString[i] == c)
                {
                    cnt++;
                    if (cnt == n)
                    {
                        result = i;
                        break;
                    }
                }
            }
            return result;
        }

    }
}
