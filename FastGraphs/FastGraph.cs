namespace FastGraphs
{
    public static class FastGraph
    {
        #region methods

        /// <summary>
        /// Read matrix from file
        /// </summary>
        /// <param name="path"> Path to file </param>
        /// <returns> Probability matrix </returns>
        static public double[,]? ReadFromFileAsync(string path)
        {
            try
            {
                double[,] matrix;
                int length = 0;
                int count = 0;
                string[] lines;
                using (StreamReader reader = new StreamReader(path))
                {
                    lines = reader.ReadToEnd().Split(new char[] { '\n' });
                    length = Int32.Parse(lines[0]);
                    count = Int32.Parse(lines[lines.Count() - 1]);
                    matrix = new double[length, length];
                    if (count != lines.Count() - 2) { throw new FileLoadException(); }
                    for (int k = 1; k < lines.Count() - 1; k++)
                    {
                        string[] line = lines[k].Split(new char[] { ' ' });
                        int i = Int32.Parse(line[0]);
                        int j = Int32.Parse(line[1]);
                        matrix[i, j] = Double.Parse(line[2]);
                        matrix[j, i] = Double.Parse(line[2]);
                    }
                    return matrix;
                }
            }
            catch
            {
                //file problems
                return null;
            }
        }

        /// <summary>
        /// Save to file
        /// </summary>
        /// <param name="path">Path to file</param>
        /// <param name="matrix">Matrix to be saved</param>
        static public void SaveFile(string path, double[,] matrix)
        {
            try
            {
                int length = matrix.GetLength(0);
                int count = 0;
                string result;
                using (StreamWriter writer = new StreamWriter(path, false))
                {
                    result = length.ToString() + '\n';
                    for (int i = 1; i < length; i++)
                    {
                        for (int j = 0; j < i; j++)
                        {
                            if (matrix[i, j] != 0)
                            {
                                result += i.ToString() + " " + j.ToString() + " " + matrix[i, j] + '\n';
                                count++;
                            }
                        }
                    }
                    result += count.ToString();
                    writer.Write(result);
                }
            }
            catch
            {
                //file problems
            }
        }

        /// <summary>
        /// Create random matrix
        /// </summary>
        /// <param name="length"></param>
        /// <returns>Matrix</returns>
        static public double[,]? GenerateRandomMatrix(int length)
        {
            try
            {
                Random random = new Random();
                double[,] matrix = new double[length, length];
                for (int i = 0; i < matrix.Length; i++)
                {
                    for (int j = i + 1; j < length; j++)
                    {
                        matrix[i, j] = random.NextDouble();
                        matrix[j, i] = matrix[i, j];
                    }
                }
                return matrix;
            }
            catch
            {
                //file problems
                return null;
            }
        }

        /// <summary>
        /// Counts degree of vertices.
        /// And also finds the maximum and minimum degree in the graph.
        /// </summary>
        static public void ExtremumDegrees(ref double[,] matrix, out int[] degreeses, out int MaxDegree, out int MaxDegreeNumber, out int MinDegree, out int MinDegreeNumber)
        {
            MinDegree = int.MaxValue;
            MaxDegree = int.MinValue;
            MaxDegreeNumber = 0;
            MinDegreeNumber = 0;
            degreeses = new int[matrix.GetLength(0)];
            for (int i = 0; i < degreeses.GetLength(0); i++)
            {
                int CurDegree = 0;
                for (int j = 0; j < degreeses.GetLength(0); j++)
                {
                    if (matrix[i, j] != 0) CurDegree++;
                }
                degreeses[i] = CurDegree;
                if (MaxDegree <= degreeses[i])
                {
                    MaxDegree = CurDegree;
                    MaxDegreeNumber = i;
                }
                if (degreeses[i] <= MinDegree)
                {
                    MinDegree = degreeses[i];
                    MinDegreeNumber = i;
                }

            }
        }

        /// <summary>
        /// Give block from matrix
        /// </summary>
        /// <param name="left">left border</param>
        /// <param name="right">right border</param>
        /// <returns></returns>
        static public double[,]? GetBlock(ref double[,] matrix, int left, int right)
        {
            double[,] output = new double[right - left + 1, right - left + 1];
            for (int i = 0; i < right - left; i++)
            {
                for (int j = 0; j < right - left; j++)
                {
                    output[i, j] = matrix[right - (right - left) + i, right - (right - left) + j];
                }
            }
            return output;
        }

        /// <summary>
        /// Finds the maximum and minimum degree in the graph. Based on degreeses.
        /// </summary>
        static public void MinMaxDegree(ref int[] degreeses, out int MaxDegree, out int MaxDegreeNumber, out int MinDegree, out int MinDegreeNumber)
        {
            MinDegree = degreeses[0];
            MaxDegree = degreeses[0];
            MaxDegreeNumber = 0;
            MinDegreeNumber = 0;
            for (int i = 1; i < degreeses.GetLength(0); i++)
            {
                if (degreeses[i] > MaxDegree)
                {
                    MaxDegree = degreeses[i];
                    MaxDegreeNumber = i;
                }
                if (degreeses[i] < MinDegree)
                {
                    MinDegree = degreeses[i];
                    MinDegreeNumber = i;
                }
            }
        }

        /// <summary>
        /// Mutual renumberin of pair of nodes.
        /// </summary>
        /// <param name="N">The first vertex of renumbering</param>
        /// <param name="M">The second vertex of renumbering</param>
        static public void Renumber(ref double[,] matrix, ref int[] degreeses, int N, int M, bool ChangeDegreses)
        {
            if (N != M)
            {
                for (int i = 0; i < matrix.GetLength(0); i++)
                {
                    double current = matrix[N, i];
                    matrix[N, i] = matrix[M, i];
                    matrix[M, i] = current;
                }
                for (int i = 0; i < matrix.GetLength(0); i++)
                {
                    double current = matrix[i, N];
                    matrix[i, N] = matrix[i, M];
                    matrix[i, M] = current;
                }
                int b = degreeses[N];
                degreeses[N] = degreeses[M];
                degreeses[M] = b;
            }
        }

        /// <summary>
        /// Checking a graph for connectivity.
        /// </summary>
        /// <returns>Returned data - connected, null - not connected.</returns>
        static public bool Сonnectivity(ref double[,] matrix)
        {
            //int[] Card = new int[matrix.GetLength(0)];
            //for (int i = 0; i < Card.GetLength(0); i++)
            //{
            //    Card[i] = i;
            //}
            //Card[0] = first;
            //Card[first] = 1;
            //int counter = 0;
            //int k = 0;
            //while ((k < matrix.GetLength(0)) & (k < counter) & (counter < matrix.GetLength(0)))
            //{
            //    for (int i = counter; i < matrix.GetLength(0); i++)
            //    {
            //        if (matrix[Card[k], Card[i]] <= 0)
            //        {
            //            counter++;
            //            int l = Card[counter];
            //            Card[counter] = Card[i];
            //            Card[i] = l;
            //        }
            //        k++;
            //    }
            //}
            //if (counter >= matrix.GetLength(0))
            //{
            //    for (int i=0;i<counter;i++)
            //    {
            //        Card[i] = 0;
            //    }
            //    return Card;
            //} 
            //else return null;





            int k = 0; int l = 0;
            int[] Mech = new int[matrix.GetLength(0)];
            Mech[0] = matrix.GetLength(0) - 1;

            bool cycle;
            while (k >= l)
            {
                for (int i = 0; i < matrix.GetLength(0); i++)
                {
                    if (matrix[Mech[l], i] > 0)
                    {
                        cycle = true;
                        int j = 0;
                        while ((j <= k) && (cycle))
                        {
                            if (i == Mech[j]) cycle = false;
                            j++;
                        }
                        if (cycle)
                        {
                            k++;
                            if (k == matrix.GetLength(0) - 1)
                            {
                                return true;
                            }
                            Mech[k] = i;
                        }
                    }
                }
                l++;
            }
            return false;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="vector"></param>
        /// <param name="ind1"></param>
        /// <param name="ind2"></param>
        /// <returns></returns>
        static public int[]? CardVectorSlice(int[] vector, int ind1, int ind2)
        {
            int[] Slice = new int[ind2 - ind1];
            if ((ind1 > ind2) || (ind2 > ind2 - ind1))
            {
                return null;
            }
            for (int i = ind1; i < ind2; i++)
            {
                Slice[i - ind1] = vector[i];
            }
            return Slice;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="vector"></param>
        /// <param name="val"></param>
        /// <returns></returns>
        static public int? InCardVector(int[] vector, int val)
        {
            if (vector.GetLength(0) == 0) return null;
            bool OK = false;
            int numb = 0;
            int i = 0;
            while (i < vector.GetLength(0) && !OK)
            {
                if (vector[i] == val)
                {
                    OK = true;
                    numb = i;
                }
                i++;
            }
            return numb;
        }
        /// <summary>
        /// Counts only degree of vertices.
        /// </summary>
        /// <param name="matrix"></param>
        /// <returns></returns>
        static public int[]? DegreeCount(double[,] matrix)
        {
            int[] degreeses = new int[matrix.GetLength(0)];
            for (int i = 0; i < degreeses.GetLength(0); i++)
            {
                int CurDegree = 0;
                for (int j = 0; j < degreeses.GetLength(0); j++)
                {
                    if (matrix[i, j] != 0) CurDegree++;
                }
                degreeses[i] = CurDegree;
            }
            return degreeses;
        }


        /// <summary>
        /// Method that pulling two degrees into one. Based on a sellected edge.
        /// </summary>
        /// <param name="mat">Probability matrix</param>
        /// <param name="i">Degree of pulling</param>
        /// <param name="j">Degree for pulling</param>
        /// <returns></returns>
        static public double[,]? EdgePull(double[,] mat, int i, int j)
        {
            if (mat[i, j] == 0) return null;
            for (int k = 0; k < mat.GetLength(0); k++)
            {
                mat[j, k] = 1 - ((1 - mat[i, k]) * (1 - mat[j, k]));
                mat[k, j] = mat[j, k];
            }
            for (int k = 0; k < mat.GetLength(0); k++)
            {
                mat[i, k] = 0;
                mat[k, i] = 0;
            }
            double[,] cur = new double[mat.GetLength(0) - 1, mat.GetLength(0) - 1];
            for (int k = 0; k < cur.GetLength(0); k++)
            {
                for (int l = 0; l < cur.GetLength(0); l++)
                {
                    cur[k, l] = mat[k, l];
                }
            }
            return cur;
        }

        ///////////

        static public void ChainReduction(ref double[,] matrix, ref int[] degreeses, int[] chain)
        {
            double Pro = 1;
            for (int i = 0; i < chain.Length; i++)
            {
                Pro *= matrix[chain[i], chain[i + 1]];
            }
            int m = matrix.Length - chain.Length;
            int[] NewNumber = new int[matrix.Length];
            NewNumber[1] = m;
            for (int i = 1; i < chain.Length; i++)
            {
                NewNumber[i] = matrix.Length - i + 1;
            }
            NewNumber[chain.Length] = matrix.Length - chain.Length;
            int NN = chain.Length;
            RenumberNodes(ref matrix, ref degreeses, chain, NewNumber, true, true);
            double[,]? matrixd = GetBlock(ref matrix, 0, m);
            int[]? degreesesd = new int[matrix.Length];
            degreesesd = CardVectorSlice(degreeses, 1, m);
            int[]? degreesesc = new int[matrix.Length];
            degreesesc = CardVectorSlice(degreeses, 1, m);
            degreesesd[m - 1]--;
            degreesesd[m]--;
            var r = matrixd[m - 1, m];
            matrixd[m - 1, m] = 0;
            matrixd[m, m - 1] = 0;
            if (r > 0)
            {
                degreesesd[m - 1]--;
                degreesesd[m]--;
            }
            int PPd = m;
            if (r > 0)
            {
                degreesesc[m - 1] = degreesesc[m] + degreesesc[m - 1] - 4;
            }
            else
            {
                degreesesc[m - 1] = degreesesc[m] + degreeses[m - 1] - 2;
            }
            int QQ = 0;
            for (int j = 0; j < m; j++)
            {
                double a = matrix[m - 1, j];
                double b = matrix[PPd, j];
                if ((a > 0) && (b > 0))
                {
                    degreesesc[m - 1]--;
                    degreesesc[j]--;
                    QQ++;
                }
                matrix[m - 1, j] = a + b - a * b;
                matrix[j, m - 1] = a + b - a * b;
            }
            matrix[m - 1, m - 1] = 0;
            var matrixc = GetBlock(ref matrix, 0, m - 1);
            int PPc = m - 1;
        }

        /// <summary>
        /// Finds a Chain if WE KNOW FOR SHURE that there is a node with degree 2
        /// </summary>
        /// <param name="matrix">Link to matrix</param>
        /// <param name="degreeses">Link to array of degreeses</param>
        /// <param name="first">First node</param>
        /// <returns>Array of chain</returns>
        static public void FindChain(ref double[,] matrix, ref int[] degreeses, out int[] Nodes, out double[] Edge_Probs, ref bool IsChain, ref bool cycle)
        {
            {
                //int newPrev = 0;

                //bool IsChain = false;
                //bool cycle = false;
                //int i = 0;

                //int[] after = new int[matrix.Length];
                //after[0] = first;
                //int prev = 1;
                //while (matrix[first, prev] == 0)
                //{ prev++; }
                //int[] befor = new int[matrix.Length];
                //befor[0] = prev;
                //int pcounter = 1;
                //int next = prev++;
                //while (matrix[first, next] == 0)
                //{ next++; }
                //after[1] = next;
                //int ncounter = 2;
                //int oldPrev = first;
                //while ((degreeses[prev] == 2) && (!cycle))
                //{
                //    newPrev = 0;
                //    while (matrix[prev, newPrev] == 0)
                //    { newPrev++; }
                //    if (newPrev == oldPrev)
                //    {
                //        newPrev++;
                //        while (matrix[prev, newPrev] == 0)
                //        { newPrev++; }
                //    }
                //    oldPrev = prev;
                //    prev = newPrev;
                //    pcounter++;
                //    befor[pcounter] = prev;
                //    if (prev == next) cycle = true;
                //}
                //int oldNext = first;
                //while (!cycle && degreeses[next] == 0)
                //{
                //    int newNext = 0;
                //    while (matrix[next, newNext] == 0)
                //    { newNext++; }
                //    if (newNext == oldNext)
                //    {
                //        newNext++;
                //        while (matrix[next, newNext] == 0)
                //        { newNext++; }
                //    }
                //    oldNext = next;
                //    next = newNext;
                //    ncounter++;
                //    after[ncounter] = next;
                //    if (prev == next) { cycle = true; }
                //}
                //int ChainLength = pcounter + ncounter;
                //int[] Nodes = new int[ChainLength];
                //for (int i = 0; i < pcounter; i++)
                //{
                //    Nodes[i] = befor[pcounter - i + 1];
                //}
                //for (int i = 0; i < ncounter; i++)
                //{
                //    Nodes[i + pcounter] = after[i];
                //}
                //double[] Probs = new double[ChainLength];
                //for (int i = 0; i < ChainLength; i++)
                //{
                //    Probs[i] = matrix[Nodes[i], Nodes[i + 1]];
                //}
                //return Probs;



                IsChain = false;
                cycle = false;
                int i = 0;

                int[] befor = null;
                int[] after = null;
                int prev, prev_count, next, next_count, old_prev, new_prev, old_next, new_next;
                while ((degreeses[i] != 2) && (i < matrix.GetLength(0) - 1))
                { i++; }
                if (i <= matrix.GetLength(0))
                {
                    befor = new int[matrix.GetLength(0) - 3];
                    after = new int[matrix.GetLength(0) - 2];
                    after[0] = i;
                    IsChain = true;
                    prev = 0;

                    while (matrix[i, prev] == 0) prev++;
                    befor[0] = prev;
                    prev_count = 1;
                    next = prev + 1;

                    while (matrix[i, next] == 0) next++;
                    after[2] = next;
                    next_count = 2;
                    //поиск предидущих

                    old_prev = 1;

                    while (degreeses[prev] == 2)
                    {
                        new_prev = 1;

                        while (matrix[prev, new_prev] == 0) new_prev++;

                        if (new_prev == old_prev)
                        {
                            old_prev = prev;
                            prev = new_prev;
                            prev_count++;
                            befor[prev_count] = prev;
                        }
                    }
                    old_next = i;

                    while (degreeses[next] == 2)
                    {
                        new_next = i;

                        while (new_next == old_next) new_next++;

                        if (new_next == old_next)
                        {
                            new_next++;
                            while (matrix[next, new_next] == 0) { new_next++; }
                        }
                        old_next = next;
                        next = new_next;
                        next_count++;
                        after[next_count] = next;
                    }
                }


                //поиск последующих закончен
                if (befor != null)
                    for (int j = 0; j < befor.GetLength(0); j++)
                    {
                        befor[j] = befor[j] + befor[befor.GetLength(0) - 1 - j];
                        befor[befor.GetLength(0) - 1 - j] = befor[j] - befor[befor.GetLength(0) - 1 - j];
                        befor[j] = befor[j] - befor[befor.GetLength(0) - 1 - j];
                    }
                Nodes = new int[befor.GetLength(0) - 1 + after.GetLength(0) - 1];
                for (int j = 0; j < befor.GetLength(0) + after.GetLength(0) - 1; j++)
                {
                    if (j < befor.GetLength(0))
                        Nodes[j] = befor[j];
                    else Nodes[j - 1] = after[j - befor.GetLength(0)];
                }
                i = 0; prev = Nodes[0];
                Edge_Probs = new double[Nodes.GetLength(0)];
                while (i < Edge_Probs.GetLength(0))
                {
                    i++;
                    next = Nodes[i - 1];
                    Edge_Probs[i - 1] = matrix[prev, next];
                    prev = next;
                }
                if (Nodes[0] == Nodes.Last()) { cycle = true; }
                else { Nodes = null; Edge_Probs = null; }
                return;
            }
        }
        static public void RenumberNodes(ref double[,] matrix, ref int[] degreeses, int[] Old, int[] New, bool Make, bool ChangeDegreeses)
        {
            //обратить внимание на исходник с N
            double[] Probs = new double[matrix.Length];
            int[] OU = new int[matrix.Length];
            int[] NewOU = new int[matrix.Length];
            int[] VIO = new int[matrix.Length];
            int[] NewVIO = new int[matrix.Length];
            var oldDegree = degreeses.DeepCopy();
            int x = 0;
            int xx = 0;

            for (int i = 0; i < matrix.Length; i++)
            {
                for (int j = 0; j < matrix.Length; j++)
                {
                    if (matrix[i, j] > 0)
                    {
                        xx++;
                        OU[xx] = i;
                        VIO[xx] = j;
                        Probs[xx] = matrix[i, j];
                    }
                }
            }
            int nCap = xx;
            if (Make)
            {
                x = 0;
                xx = Old.Length;
                for (int i = 0; i < matrix.Length; i++)
                {
                    if ((InCardVector(Old, i) == null) && (InCardVector(New, i) != null))
                    {
                        xx++;
                        Old[xx] = i;
                        while ((InCardVector(Old, x) == null) && (InCardVector(New, x) != null)) { x++; }
                        New[xx] = x;
                        x++;
                    }
                }
                //N = xx;
            }
            if (ChangeDegreeses)
            {
                for (int i = 0; i < degreeses.GetLength(0); i++)
                {
                    degreeses[New[i]] = oldDegree[Old[i]];
                }
            }
        }

        static public double Reduction(ref double[,] mat, ref int[] degreeses, int length)
        {
            //double R = 0, dL, dU = 0, WW, RP=0;
            //int m=0, iii;
            //int minDegree;
            //int minDegreeNumber;
            //int maxDegree;
            //int maxDegreeNumber;
            //MinMaxDegree(ref degreeses,out maxDegree, out maxDegreeNumber, out minDegree, out minDegreeNumber);
            //while ((minDegree==1) && (matrix.Length>2))
            //{
            //    if (minDegreeNumber<=matrix.Length)
            //    {
            //        m = 1;
            //        while ((m<=matrix.Length) && (matrix[minDegreeNumber, m]==0))
            //        {
            //            m++;
            //        }
            //        RP = matrix[minDegreeNumber, m];
            //        degreeses[m]--;
            //        iii = degreeses[matrix.Length];
            //        Renumber(ref matrix, ref degreeses, matrix.Length,minDegreeNumber, true);
            //        degreeses[minDegreeNumber] = iii;
            //    }
            //    ;
            //    length--;
            //    MinMaxDegree(ref degreeses, out maxDegree, out maxDegreeNumber, out minDegree, out minDegreeNumber);
            //    dL = Pr * R;
            //    dU = dU * Pr;
            //    Fin = true;



            int ii;
            bool Fin = false;
            if (!FastGraph.Сonnectivity(ref mat))
            { return 0; }
            // Graph optimization
            int MaxDegree = 0; int MaxDegreeNumber = 0; int MinDegreeNumber = 0; int minDeg = 1; double R = 1; int m = 0; double RP = 0;
            // dangling vertex removal
            while ((minDeg == 1) && (mat.GetLength(0) > 1))
            {
                FastGraph.MinMaxDegree(ref degreeses, out MaxDegree, out MaxDegreeNumber, out minDeg, out MinDegreeNumber);
                if ((minDeg == 1) && (MinDegreeNumber <= mat.GetLength(0)))
                {
                    for (int i = 0; i < mat.GetLength(0); i++)
                    {
                        if (mat[MinDegreeNumber, i] > 0)
                        { m = i; } //number of adjacent vertex
                        RP = mat[MinDegreeNumber, m];
                        ii = degreeses[m];
                        FastGraph.Renumber(ref mat, ref degreeses, MinDegreeNumber, mat.Length - 1, true);
                        R = R * RP;
                        mat = FastGraph.GetBlock(ref mat, 0, mat.GetLength(0) - 1);
                    }
                }
            }

            //dimension 1

            if (mat.GetLength(0) == 1)
            { Fin = true; return R; }


            //cycle
            if (minDeg == 2)
            {
                if ((minDeg == MaxDegree) && (minDeg == 2))
                {
                    List<double> list = new List<double>();
                    double RR = 1;
                    ii = 1; int x = 0;
                    for (int i = 0; i < mat.GetLength(0); i++)
                    {
                        int iii = 1;
                        while ((mat[ii, iii] == 0) && (iii == x))
                        { iii++; }
                        list[i] = (mat[ii, iii]);
                        RR = RR * list[i];
                        x = ii;
                        ii = iii;
                    }
                    RP = (1 - list[0]) / list[0];
                    for (int i = 1; i < mat.GetLength(0); i++)
                    {
                        RP = RP + ((1 - list[i]) / list[i]);
                    }
                    RP = RR * (1 + RP);

                    R = R * RP; Fin = true; return R;
                }

                int ChainCounter = 0;
                double CL = 0;
                bool IsCycle = false;
                bool IsChain = true;
                double[] ch;
                int[] Nodes;
                while ((IsChain) && mat.GetLength(0) > 4)
                {
                    FastGraph.FindChain(ref mat, ref degreeses, out Nodes, out ch, ref IsChain, ref IsCycle);
                    if (Nodes != null)
                    {
                        ChainCounter++;
                        CL += Nodes.GetLength(0);
                        if (!IsCycle)
                        {
                            m = mat.GetLength(0) - Nodes.GetLength(0);
                            int[] NewNumbers = new int[Nodes.GetLength(0)];
                            NewNumbers[0] = mat.GetLength(0) - Nodes.GetLength(0);
                            for (int i = 1; i < Nodes.GetLength(0) - 1; i++)
                            { NewNumbers[i] = mat.GetLength(0) - i + 1; }
                            NewNumbers[Nodes.GetLength(0)] = mat.GetLength(0) - Nodes.GetLength(0);
                            FastGraph.RenumberNodes(ref mat, ref degreeses, Nodes, NewNumbers, true, true);
                            RP = 0;
                            for (int i = 0; i < Nodes.GetLength(0); i++)
                            { RP += 1 / ch[i]; }
                            RP = 1 / (RP - Nodes.GetLength(0));
                            double LastEdgeProb = ch[Nodes.GetLength(0)];

                            mat = FastGraph.GetBlock(ref mat, 0, m - 1);

                            double r = mat[m - 1, m];
                            mat[m - 1, m] = r + RP - r * RP;
                            mat[m, m - 1] = r + RP - r * RP;

                            //if (r > 0) { degreeses}
                        }
                        else
                        {
                            double[] CH = new double[Nodes.GetLength(0)];
                            for (int i = 0; i < length - 2; i++) { CH[i] = ch[i]; }
                            ch = CH;
                            double RR = 1;
                            for (int i = 0; i < Nodes.GetLength(0); i++) { RR *= ch[i]; }
                            RP = (1 - ch[1]) / ch[1];
                            for (int i = 1; i < ch.GetLength(0); i++) { RP = RP + (1 - ch[i]) / ch[i]; }
                            RP = RR * (1 + RP);
                            int[] NewNumbers = new int[ch.GetLength(0)];
                            NewNumbers[0] = mat.GetLength(0) - ch.GetLength(0);
                            for (int i = 1; i < ch.GetLength(0) - 1; i++)
                            { NewNumbers[i] = mat.GetLength(0) - i + 1; }
                            mat = FastGraph.GetBlock(ref mat, 0, mat.GetLength(0) - ch.GetLength(0) + 1);
                            R = R * RP;
                        }
                    }
                    else IsChain = false;
                }
            }
            return R;
        }
        //////////

        #endregion
    }
}