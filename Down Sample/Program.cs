using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DownSample
{
    class Program
    {
        public static utl.GenMapV3 gMap;
        static void Main(string[] args)
        {
            string gMap_Path = "";
            string seq_VCF_Path = "";
            string seqVCF_statFilePath = "";
            string arrVCF_collectionList = "";
            string seqVCF_OutputPath = "";

            int target_nSite = 17100;
            int windowSize = 5;


            gMap = new utl.GenMapV3(gMap_Path,1,3);

            //Step 1 is a long process, we only need to run once.
            //There is parallel option within to adjuct for more cpu cores
            DownSample.Step1 sp1 = new DownSample.Step1();
            sp1.Run(seq_VCF_Path, seqVCF_statFilePath);

            //Step 2 generate a collection by using a window size.
            //This module will output a # of marker to be collected by using this windows size.
            //Try different windows sizes to reach desired # of marker.
            DownSample.Step2 sp2 = new DownSample.Step2();
            sp2.run(target_nSite, windowSize, seqVCF_statFilePath, arrVCF_collectionList);

            //Step 3 will simply extract array data from the sequencing data by list of site generated from step 2.
            DownSample.Step3 sp3 = new DownSample.Step3();
            sp3.Run(arrVCF_collectionList, seqVCF_statFilePath, seqVCF_OutputPath);


        }
    }
}
