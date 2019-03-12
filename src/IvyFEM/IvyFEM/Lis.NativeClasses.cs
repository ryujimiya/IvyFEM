using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using LisScalar = System.Double;

namespace IvyFEM.Lis
{
    unsafe public struct NativeLisCommTable
    {
        private int comm;
        private int pad;
        private int neibpetot;
        private int imnnz;
        private int exnnz;
        private int wssize;
        private int wrsize;
        private int* neibpe;
        private int* import_ptr;
        private int* import_index;
        private int* export_ptr;
        private int* export_index;
        private LisScalar* ws;
        private LisScalar* wr;
        //private MPI_Request* req1;
        //private MPI_Request* req2;
        //private MPI_Status* sta1;
        //private MPI_Status* sta2;
    }

    unsafe public struct NativeLisVector
    {
        private int label;
        private int status;
        private int precision;
        private int gn;
        private int n;
        private int np;
        private int pad;
        private int origin;
        private int is_copy;
        private int is_destroy;
        private int is_scaled;
        private int my_rank;
        private int nprocs;
        private int comm;
        private int @is;
        private int ie;
        private int* ranges;
        private LisScalar* value;
        private LisScalar* value_lo;
        private LisScalar* work;
        private int intvalue;
    }

    unsafe public struct NativeLisMatrixCore
    {
        private int nnz;
        private int ndz;
        private int bnr;
        private int bnc;
        private int nr;
        private int nc;
        private int bnnz;
        private int nnd;
        private int maxnzr;
        private int* ptr;
        private int* row;
        private int* col;
        private int* index;
        private int* bptr;
        private int* bindex;
        private LisScalar* value;
        private LisScalar* work;
    }

    unsafe public struct NativeLisMatrixDiag
    {
        private int label;
        private int status;
        private int precision;
        private int gn;
        private int n;
        private int np;
        private int pad;
        private int origin;
        private int is_copy;
        private int is_destroy;
        private int is_scaled;
        private int my_rank;
        private int nprocs;
        private int comm;
        private int @is;
        private int ie;
        private int* ranges;
        private LisScalar* value;
        private LisScalar* work;

        private int bn;
        private int nr;
        private int* bns;
        private int* ptr;
        private LisScalar** v_value;
    }

    unsafe public struct NativeLisMatrix
    {
        private int label;
        private int status;
        private int precision;
        private int gn;
        private int n;
        private int np;
        private int pad;
        private int origin;
        private int is_copy;
        private int is_destroy;
        private int is_scaled;
        private int my_rank;
        private int nprocs;
        private int comm;
        private int @is;
        private int ie;
        private int* ranges;

        private int matrix_type;
        private int nnz; /* CSR,CSC,MSR,JAD,VBR,COO */
        private int ndz; /* MSR */
        private int bnr; /* BSR,BSC */
        private int bnc; /* BSR,BSC */
        private int nr; /* BSR,BSC,VBR */
        private int nc; /* BSR,BSC,VBR */
        private int bnnz; /* BSR,BSC,VBR */
        private int nnd; /* DIA */
        private int maxnzr; /* ELL,JAD */
        private int* ptr; /* CSR,CSC,JAD */
        private int* row; /* JAD,VBR,COO */
        private int* col; /* JAD,VBR,COO */
        private int* index; /* CSR,CSC,MSR,DIA,ELL,JAD */
        private int* bptr; /* BSR,BSC,VBR */
        private int* bindex; /* BSR,BSC,VBR */
        private LisScalar* value; /* CSR,CSC,MSR,DIA,ELL,JAD,BSR,BSC,VBR,DNS,COO */
        private LisScalar* work;

        private NativeLisMatrixCore* L;
        private NativeLisMatrixCore* U;
        private NativeLisMatrixDiag* D;
        private NativeLisMatrixDiag* WD;

        private int is_block;
        private int pad_comm;
        private int is_pmat;
        private int is_sorted;
        private int is_splited;
        private int is_save;
        private int is_comm;
        private int is_fallocated;
        private int use_wd;
        private int conv_bnr;
        private int conv_bnc;
        private int* conv_row;
        private int* conv_col;
        //int[] options = new int[IvyFEM.Lis.Constants.LisMatrixOptionLen];
        private int option1;
        private int option2;
        private int option3;
        private int option4;
        private int option5;
        private int option6;
        private int option7;
        private int option8;
        private int option9;
        private int option10;

        private int w_annz;
        private int* w_nnz;
        private int* w_row;
        private int** w_index;
        private LisScalar** w_value;
        private LisScalar*** v_value;

        private int* l2g_map;
        private NativeLisCommTable* commtable;
    }

    unsafe public struct NativeLisMatrixILU
    {
        private int n;
        private int bs;
        private int* nnz_ma;
        private int* nnz;
        private int* bsz;
        private int** index;
        private LisScalar** value;
        private LisScalar*** values;
    }


    unsafe public struct NativeLisPrecon
    {
        private int precon_type;
        private NativeLisMatrix* A; /* SSOR */
        private NativeLisMatrix* Ah;
        private NativeLisMatrixILU* L; /* ilu(k),ilut,iluc,sainv */
        private NativeLisMatrixILU* U; /* ilu(k),ilut,iluc,sainv */
        private NativeLisMatrixDiag* WD; /* bilu(k),bilut,biluc,bjacobi */
        private NativeLisVector* D; /* ilu(k),ilut,iluc,jacobi,sainv */
        private NativeLisVector* Pb; /* i+s */
        private NativeLisVector* temp; /* saamg */
        private double theta; /* saamg */
        private NativeLisVector** work; /* adds */
        private NativeLisSolver* solver; /* hybrid */ // LisSolver *
        private int worklen; /* adds */
        private int level_num; /* saamg */
        private int wsize; /* saamg */
        private int solver_comm; /* saamg */
        private int my_rank; /* saamg */
        private int nprocs; /* saamg */
        private int is_copy;
        private NativeLisCommTable* commtable; /* saamg */
    }

    unsafe public struct NativeLisSolver
    {
        private NativeLisMatrix* A;
        private NativeLisMatrix* Ah;
        private NativeLisVector* b;
        private NativeLisVector* x;
        private NativeLisVector* xx;
        private NativeLisVector* d;
        private NativeLisMatrixDiag* WD;
        private NativeLisPrecon* precon;
        private NativeLisVector** work;
        private double* rhistory;
        private int worklen;
        //int[] options = new int[IvyFEM.Lis.Constants.LisOptionsLen];
        private int option1;
        private int option2;
        private int option3;
        private int option4;
        private int option5;
        private int option6;
        private int option7;
        private int option8;
        private int option9;
        private int option10;
        private int option11;
        private int option12;
        private int option13;
        private int option14;
        private int option15;
        private int option16;
        private int option17;
        private int option18;
        private int option19;
        private int option20;
        private int option21;
        private int option22;
        private int option23;
        private int option24;
        private int option25;
        private int option26;
        private int option27;
        //LisScalar[] @params = new double[IvyFEM.Lis.Constants.LisParamsLen];
        private LisScalar param1;
        private LisScalar param2;
        private LisScalar param3;
        private LisScalar param4;
        private LisScalar param5;
        private LisScalar param6;
        private LisScalar param7;
        private LisScalar param8;
        private LisScalar param9;
        private LisScalar param10;
        private LisScalar param11;
        private LisScalar param12;
        private LisScalar param13;
        private LisScalar param14;
        private LisScalar param15;
        private int retcode;
        private int iter;
        private int iter2;
        private double resid;
        private double time;
        private double itime;
        private double ptime;
        private double p_c_time;
        private double p_i_time;
        private int precision;
        private double bnrm;
        private double tol;
        private double tol_switch;
        private int setup;
    }
}
