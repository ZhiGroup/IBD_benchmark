using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Collections.Concurrent;


namespace DownSample
{
    class DownSample:Program
    {
        /// <summary>
        /// run a phy gen minor count scan first
        /// use a VCF input then output a file contains phyical location, genetic location, and minor allele count
        /// </summary>
        public class Step1
        {

            public static int nWorkerThread = 10;//parallel run,put more number if you have more cpu core.
            public static int scale = 10000;//progress tracking scale.
            public static int msWorkerPoolCheck = 5000;
            public static int msReaderWait = 100;
            public static int msWorkerWait = 100;
            ConcurrentDictionary<int, bool> WorkerPool = new ConcurrentDictionary<int, bool>();
            BlockingCollection<dataRow> procQueue;
            public static Phy_Gen_MiAC[] holder;
            public static int nWorker = 10;
     
            /// <summary>
            /// output phy, gen, MiAC
            /// gives minor allele count -1 if the row is Multi allele.
            /// 
            /// If the nSiteInFile is not given (if we don't now) the program will run a scan to know how many sites in the file.
            /// Simulator usually output # of sites created in the log.
            /// </summary>
            /// <param name="inPath">sequencing VCF path</param>
            /// <param name="outPath"></param>
            /// <param name="nSiteInFile"># of site in sequencing data</param>
            /// <param name="buffSize"></param>
            public void Run(string inPath, string outPath, int nSiteInFile=0, int buffSize = 50000)
            {
                procQueue = new BlockingCollection<dataRow>(buffSize);
   
                nWorker = nWorkerThread;
                Console.WriteLine(DateTime.Now + " Making Holder...");

                if (nSiteInFile == 0)
                {
                    nSiteInFile = utl.lineCount(inPath, "#");
                }

                holder = new Phy_Gen_MiAC[nSiteInFile];
                Console.WriteLine(DateTime.Now + " Holder Created.");

                DateTime st = DateTime.Now;
                runReader(inPath, nSiteInFile);

                runWorkers();


                Wait();

                TimeSpan ts = DateTime.Now - st;
                Console.WriteLine(ts.TotalMilliseconds);


                WriteResult(outPath);

            }


            void WriteResult(string outPath)
            {
                Console.WriteLine(DateTime.Now + "Writing Results...");
                StreamWriter sw = new StreamWriter(outPath);
                sw.NewLine = "\n";
                foreach (Phy_Gen_MiAC one in holder)
                {
                    sw.WriteLine(one.phy + "\t" + one.gen + "\t" + one.miAC);
                }


                sw.Close();
            }

            void Reader(string inPath, int nSiteInFile)
            {

                StreamReader sr = new StreamReader(inPath);
                string line;

                while ((line = sr.ReadLine()) != null && line.StartsWith("##"))
                {
                    continue;
                }
                Console.WriteLine(DateTime.Now + " " + inPath + " Reader Started.");
                int rowCnt = 0;
                /////////////////////////////////////////////////////////

                DateTime sTimeOut = DateTime.Now;
                while ((line = sr.ReadLine()) != null)
                {

                    while (!procQueue.TryAdd(new dataRow(rowCnt, line)))
                    {
                        System.Threading.Thread.Sleep(msReaderWait);
                    }

                    rowCnt++;

                }
                sr.Close();
                TimeSpan ts = DateTime.Now - sTimeOut;
                //Console.WriteLine("Queue full " + procQueue.Count() + "\t" + ts.TotalMilliseconds);
                procQueue.CompleteAdding();
                Console.WriteLine(rowCnt);

                Console.WriteLine(DateTime.Now + " " + inPath + " Reader: Read Completed.");


            }


            void runReader(string inPath, int nSiteInFile)
            {
                Task one = Task.Factory.StartNew(() => Reader(inPath, nSiteInFile));
            }


            void runWorkers()
            {
                for (int i = 0; i < nWorker; i++)
                {
                    Task one = Task.Factory.StartNew(() => Worker());
                }

            }

            void Worker()
            {
                Console.WriteLine(DateTime.Now + " Worker " + System.Threading.Thread.CurrentThread.ManagedThreadId + " Online.");
                WorkerPool.TryAdd(System.Threading.Thread.CurrentThread.ManagedThreadId, false);

                dataRow data = new dataRow(0, "");

                //int abc = procQueue.Count();
                //Console.WriteLine(abc);
                while (procQueue.Count() != 0 || procQueue.IsAddingCompleted == false)
                {

                    if (procQueue.TryTake(out data) == false)
                    {
                        System.Threading.Thread.Sleep(msWorkerWait);
                        continue;
                    }
                    string line = data.Content;
                    int rowIndex = data.RowIndexInFile;
                    holder[rowIndex] = new Phy_Gen_MiAC(line);

                }




                WorkerPool[System.Threading.Thread.CurrentThread.ManagedThreadId] = true;
                Console.WriteLine(DateTime.Now + " Worker " + System.Threading.Thread.CurrentThread.ManagedThreadId + " Offline.");

            }

            void Wait()
            {
                Console.WriteLine(DateTime.Now + " Waiting Workers...");
                bool waitComplete = false;
                while (waitComplete == false)
                {
                    System.Threading.Thread.Sleep(msWorkerPoolCheck);
                    waitComplete = true;

                    if (procQueue.IsAddingCompleted == false)
                    {
                        waitComplete = false;
                        continue;
                    }

                    foreach (bool done in WorkerPool.Values)
                    {
                        if (done == false)
                        {
                            waitComplete = false;
                            break;
                        }
                    }
                }

                Console.WriteLine(DateTime.Now + " Workers All Finish.");
            }

            class dataRow
            {
                public int RowIndexInFile;
                public string Content;
                public dataRow(int rowIndex, string line)
                {
                    RowIndexInFile = rowIndex;
                    Content = line;
                }

            }


        }

        public class Phy_Gen_MiAC
        {
            public int phy = 0;
            public double gen = -1;
            public int miAC = -1;
            public Phy_Gen_MiAC(string line)
            {
                //#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NWD112649	NWD278543
                string[] parts = line.Split('\t');
                phy = Convert.ToInt32(parts[1]);
                if (parts[3].Length != 1 || parts[4].Length != 1)
                { return; }
                gen = gMap.getGenLoc(phy);
                miAC = utl.MinorAlleleCnt(parts);
            }

            public Phy_Gen_MiAC(string line, bool Load_PGM)
            {
                string[] parts = line.Split('\t');
                phy = Convert.ToInt32(parts[0]);
                gen = Convert.ToDouble(parts[1]);
                miAC = Convert.ToInt32(parts[2]);
            }


        }
        /// <summary>
        /// down sample by the result of step1
        /// 1 load pgm file(result of step1)
        /// 2 determing which site to take create a collection list, row index of site.
        /// 3 output collection list to a file
        /// </summary>
        public class Step2
        {
            public static int hash_Fail = 0;
            public static int nWnd_PerGP = 0;
            public static int failCnt = 0;
            public static int fail_Short = 0;

            public static double sWndSize = 0;
            public static double bWndSize = 0;


            public class Index_PGM
            {
                public int index = 0;
                public Phy_Gen_MiAC PGM;
                public Index_PGM(int rowIndex, Phy_Gen_MiAC pgmData)
                {
                    index = rowIndex;
                    PGM = pgmData;
                }


            }

            public static List<Index_PGM> PGM_Holder = new List<Index_PGM>();

            /// <summary>
            /// Try different # of intervial to reach your desired # of site in new array data.
            /// Note, using bigger # of intervial (nWnd) results more uneven the markers distributed, but more markers will be selected.
            /// On the other hand, using small window will give more even markers, but less markers.
            /// 
            /// We believed using 90% of the target # of marker is good enough.
            /// </summary>
            /// <param name="nSiteArr">Target # of site in array data, we may use a reference number such as # of site in UKB array data</param>
            /// <param name="nWnd"> # of intervial window to try</param>
            /// <param name="pgmFilePath"></param>
            /// <param name="outPath"></param>
            /// <returns></returns>
            public int run(double nSiteArr, int nWnd, string pgmFilePath, string outPath)
            {
                nWnd_PerGP = nWnd;

                int gpCnt = 0;

                List<Index_PGM> holder = Load_PGM(pgmFilePath);
                sWndSize = (holder.Last().PGM.gen - holder.First().PGM.gen) / nSiteArr;
                bWndSize = sWndSize * nWnd;


                HashSet<int> rowIndexToTake = new HashSet<int>();
                #region make index hash

                List<Index_PGM> oneGP_Holder = new List<Index_PGM>();
                //List<List<Index_PGM>> gpHolder = new List<List<Index_PGM>>();
                double L_Range_Gen = holder.First().PGM.gen;
                for (int i = 0; i < holder.Count(); i++)
                {

                    if (holder[i].PGM.gen - L_Range_Gen >= bWndSize)
                    {//push a window
                        L_Range_Gen = holder[i].PGM.gen;

                        ProcessOne_GigWnd(oneGP_Holder, rowIndexToTake);
                        gpCnt++;

                        oneGP_Holder.Clear();

                        //next window
                        oneGP_Holder.Add(holder[i]);
                    }
                    else
                    {
                        oneGP_Holder.Add(holder[i]);
                    }


                }

                //handle tail
                ProcessOne_GigWnd(oneGP_Holder, rowIndexToTake);

                gpCnt++;
                #endregion
                Console.WriteLine("wSize= " + nWnd + ". Collected " + rowIndexToTake.Count() + " , Missed " + (nSiteArr - rowIndexToTake.Count()));
                //Console.WriteLine("GP Cnt " + gpCnt +" by "+nWnd_PerGP);

                utl.listToFile(rowIndexToTake.ToList(), outPath);

                return rowIndexToTake.Count();
            }

            public static void ProcessOne_GigWnd(List<Index_PGM> bigGpHolder, HashSet<int> rowIndexToTake)
            {
                List<List<Index_PGM>> gpHolder = Break_BigWnd(bigGpHolder);

                List<int> nTakes = TakeN_From_N_Wnd(gpHolder);

                for (int t = 0; t < nTakes.Count(); t++)
                {
                    List<int> real_Indexs = TakeN_From_OneWnd(nTakes[t], gpHolder[t]);
                    foreach (int one in real_Indexs)
                    {
                        rowIndexToTake.Add(one);
                    }
                }

            }

            /// <summary>
            /// break a big window into small windows
            /// </summary>
            /// <param name="bigGpHolder"></param>
            /// <returns></returns>
            public static List<List<Index_PGM>> Break_BigWnd(List<Index_PGM> bigGpHolder)
            {
                List<Index_PGM> oneGP_Holder = new List<Index_PGM>();
                List<List<Index_PGM>> gpHolder = new List<List<Index_PGM>>();
                double L_Range_Gen = bigGpHolder.First().PGM.gen;
                for (int i = 0; i < bigGpHolder.Count(); i++)
                {

                    if (bigGpHolder[i].PGM.gen - L_Range_Gen >= sWndSize)
                    {//push a window
                        L_Range_Gen = bigGpHolder[i].PGM.gen;

                        gpHolder.Add(new List<Index_PGM>(oneGP_Holder));

                        oneGP_Holder.Clear();
                        //next window
                        oneGP_Holder.Add(bigGpHolder[i]);
                    }
                    else
                    {
                        oneGP_Holder.Add(bigGpHolder[i]);
                    }

                }

                if (oneGP_Holder.Count() != 0)
                {
                    gpHolder.Add(new List<Index_PGM>(oneGP_Holder));
                }

                return gpHolder;
            }


            public static void ProcessOneGP(List<List<Index_PGM>> gpHolder, HashSet<int> rowIndexToTake)
            {

                List<int> nTakes = TakeN_From_N_Wnd(gpHolder);
                //Console.WriteLine("nTakes.Sum(): " + nTakes.Sum());
                if (nTakes.Sum() != nWnd_PerGP)
                {
                    Console.WriteLine("TakeN_From_N_Wnd Failed " + nTakes.Sum());
                    failCnt++;
                    fail_Short += (nWnd_PerGP - nTakes.Sum());
                }
                for (int t = 0; t < nTakes.Count(); t++)
                {
                    List<int> real_Indexs = TakeN_From_OneWnd(nTakes[t], gpHolder[t]);
                    foreach (int one in real_Indexs)
                    {
                        rowIndexToTake.Add(one);
                    }
                }

            }



            public static List<Index_PGM> Load_PGM(string pgmPath)
            {
                List<Index_PGM> res = new List<Index_PGM>();
                string line;

                int rowCnt = 0;
                StreamReader sr = new StreamReader(pgmPath);
                while ((line = sr.ReadLine()) != null)
                {
                    Phy_Gen_MiAC pgm = new Phy_Gen_MiAC(line, true);
                    res.Add(new Index_PGM(rowCnt, pgm));

                    rowCnt++;
                }
                sr.Close();
                return res;
            }

            /// <summary>
            /// give the final row index to take
            /// </summary>
            /// <param name="n"></param>
            /// <param name="rows"></param>
            /// <returns></returns>
            public static List<int> TakeN_From_OneWnd(int n, List<Index_PGM> rows)
            {
                List<int> res = new List<int>();
                List<List<Index_PGM>> divided_PGM = utl.DivideList(rows, n);
                foreach (List<Index_PGM> oneGP in divided_PGM)
                {
                    res.Add(TakeOne_From_OneWnd(oneGP));
                }

                return res;
            }

            /// <summary>
            /// give the final row index
            /// </summary>
            /// <param name="rows"></param>
            /// <returns></returns>
            public static int TakeOne_From_OneWnd(List<Index_PGM> rows)
            {
                int maxMiC = 0;
                int bestIndex = 0;
                foreach (Index_PGM one in rows)
                {
                    if (one.PGM.miAC > maxMiC)
                    {
                        maxMiC = one.PGM.miAC;
                        bestIndex = one.index;
                    }
                }

                return bestIndex;
            }

            /// <summary>
            /// give a list of count 
            /// </summary>
            /// <param name="counts"></param>
            /// <returns></returns>
            public static List<int> TakeN_From_N_Wnd(List<List<Index_PGM>> gpHolder)
            {
                List<int> counts = new List<int>();
                foreach (List<Index_PGM> oneGP in gpHolder)
                {
                    counts.Add(oneGP.Count());
                }



                List<int> res = new List<int>();
                for (int i = 0; i < counts.Count(); i++)
                {
                    res.Add(0);
                }

                int nTaken = 0;
                while (true)
                {
                    for (int i = 0; i < counts.Count(); i++)
                    {
                        if (counts[i] > 0)
                        {
                            res[i]++;
                            counts[i]--;
                            nTaken++;
                        }

                        if (nTaken == nWnd_PerGP)
                        {
                            return res;
                        }
                    }
                    //nothing to take break
                    if (counts.Sum() == 0)
                    {
                        return res;
                    }
                }


            }

            public void Load_PGMs(string inPath)
            {
                StreamReader sr = new StreamReader(inPath);
                string line;

                int rowIndex = 0;
                Phy_Gen_MiAC temPgm;
                while ((line = sr.ReadLine()) != null)
                {
                    temPgm = new Phy_Gen_MiAC(line);
                    if (temPgm.miAC != -1)
                    {
                        PGM_Holder.Add(new Index_PGM(rowIndex, temPgm));
                    }

                    rowIndex++;
                }
                sr.Close();
            }

        }

        /// <summary>
        /// Extract array data from sequencing data by given a list of site to take.
        /// The list of site file is generated by step2
        /// </summary>
        public class Step3
        {
            public static int scale=100000;

            void run_Inner(string listOfSite, string vcfPath, string outPath, int scale)
            {

                Console.WriteLine("Extract Array from Sequencing.");
                Console.WriteLine("listOfSite: " + listOfSite);
                Console.WriteLine("vcfPath: " + vcfPath);
                Console.WriteLine("outPath: " + outPath);
                Console.WriteLine("scale: " + scale);


                List<String> siteIndexs = utl.fileToList(listOfSite);
                HashSet<int> rowIndexToTake = new HashSet<int>(siteIndexs.Select(int.Parse).ToList());

                string line;
                StreamReader sr = new StreamReader(vcfPath);
                StreamWriter sw = new StreamWriter(outPath);
                sw.NewLine = "\n";
                while ((line = sr.ReadLine()) != null && line.StartsWith("##"))
                {
                    sw.WriteLine(line);
                }
                sw.WriteLine(line);
                int rowCnt = 0;
                while ((line = sr.ReadLine()) != null)
                {
                    if (rowIndexToTake.Contains(rowCnt))
                    {
                        sw.WriteLine(line);
                    }


                    if (rowCnt % scale == 0)
                    {
                        Console.WriteLine(DateTime.Now + " " + rowCnt);

                    }

                    rowCnt++;
                }


                sr.Close();
                sw.Close();



            }

            public void Run(string[] args)
            {
                string listOfSite = args[0];
                string vcfPath = args[1];
                string outPath = args[2];


                run_Inner(listOfSite, vcfPath, outPath, scale);

            }

            public void Run(string listOfSite, string vcfPath, string outPath)
            {
                run_Inner(listOfSite, vcfPath, outPath, scale);

            }

        }
    }
}
