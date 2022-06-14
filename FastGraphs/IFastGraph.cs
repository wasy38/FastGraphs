using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FastGraphs
{
    internal interface IFastGraph
    {
        /// <summary>
        /// Read matrix from file
        /// </summary>
        /// <param name="path"> Path to file </param>
        /// <returns> Probability matrix </returns>
        public double[,]? ReadFromFileAsync(string path);

        /// <summary>
        /// Save to file
        /// </summary>
        /// <param name="path">Path to file</param>
        /// <param name="matrix">Matrix to be saved</param>
        public void SaveFile(string path, double[,] matrix);

        /// <summary>
        /// Create random matrix
        /// </summary>
        /// <param name="length"></param>
        /// <returns>Matrix</returns>
        public double[,]? GenerateRandomMatrix(int length);

        /// <summary>
        /// Mutual renumberin of pair of nodes.
        /// </summary>
        /// <param name="N">The first vertex of renumbering</param>
        /// <param name="M">The second vertex of renumbering</param>
        public void Renumber(ref double[,] matrix, ref int[] degreeses, int N, int M, bool ChangeDegreses);

        /// <summary>
        /// Finds the maximum and minimum degree in the graph.
        /// </summary>
        public void MinMaxDegree(ref int[] degreeses, out int MaxDegree, out int MaxDegreeNumber, out int MinDegree, out int MinDegreeNumber);

        /// <summary>
        /// Give block from matrix
        /// </summary>
        /// <param name="left">left border</param>
        /// <param name="right">right border</param>
        /// <returns></returns>
        public double[,]? GetBlock(ref double[,] matrix, int left, int right);

        /// <summary>
        /// Counts degree of vertices.
        /// And also finds the maximum and minimum degree in the graph.
        /// </summary>
        public void ExtremumDegrees(ref double[,] matrix, out int[] degreeses, out int MaxDegree, out int MaxDegreeNumber, out int MinDegree, out int MinDegreeNumber);

        /// <summary>
        /// Checking a graph for connectivity.
        /// </summary>
        /// <returns>True - connected, Fale - not connected.</returns>
        public int[]? Сonnectivity(ref double[,] matrix, ref int[] degreeses, int first);

        public int? InCardVector(int[] vector, int val);

        public int[]? CardVectorSlice(int[] vector, int ind1, int ind2);

        /////////////////////

        /// <summary>
        /// Graph chain reduction
        /// </summary>
        public void ChainReduction(ref double[,] matrix, ref int[] degreeses, int[] chain);

        /// <summary>
        /// Finds a Chain if WE KNOW FOR SHURE that there is a node with degree 2
        /// </summary>
        public void FindChain(ref double[,] matrix, ref int[] degreeses, int first);

        /// <summary>
        /// Renumbering nodes
        /// </summary>
        public void RenumberNodes(ref double[,] matrix, ref int[] degreeses, int[] Old, int[] New, bool Make, bool ChangeDegreeses);

        public double Reduction(ref double[,] mat, ref int[] degreeses);
    }
}
