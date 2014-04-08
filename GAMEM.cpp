#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define MAXL 2000
#define MAXADJ 2000
#define MAXN 2000
#define MAXSLOT 35
#define MAXCOST 10000000
#define N_loop 200     //so the he
#define N_size 200     //kich thuoc quan the
#define Pc 0.2
#define Pm 0.05
#define Ps 0.3
#define alpha 2
#define eps 0.0000000001
#define MAXDELAY 10000000

char *FILEINP = "", *FILEOUT = "";
struct List
{
    long N;
    long element[MAXL];
};

struct Tree
{
    long M;
    long A[MAXN], B[MAXN];
};

struct Sche
{
    long N;
    long Node[MAXN];
    long Slot[MAXN];
};

FILE *fi, *fo;
clock_t start_time, finish_time;
long Es = 100, Er = 15;

double run_time;
long N, M, K, S;
int A[MAXN][MAXN];
int A_R[MAXN][MAXN];    //canh cua cay result

List adj[MAXN], grama[MAXN], ListTer;
List child[MAXN];
List send[MAXN];

long Ter[MAXN];
int isTer[MAXN];

Tree T[N_size + 5];
Sche B[N_size + 5];
Tree T_KQ;
Sche B_KQ;

int X_conflict[MAXN][MAXN];
long n_color_slot[MAXSLOT], color_slot[MAXSLOT][MAXN];

void process();
void input();
void khoitao();
long min_hit_set(List [], List &, List &);
void dij(List [], List [], List &, List []);
double rand01();
void ganList(List &, List &);
void ganTree(Tree &, Tree &);
void ganSche(Sche &, Sche &);
void LGRand(Tree &, Sche &, Tree &, Sche &, Tree &, Sche &);
void LGWeight(Tree &, Sche &, Tree &, Sche &, Tree &, Sche &);
void dotbien(Tree &, Sche &);
int inList(long, List &);

int main(int Nts, char **arg)
//int main()
{

    srand(time(NULL));
    FILEINP = arg[1];
    FILEOUT = arg[2];
    process();
    fclose(fi);
    fclose(fo);
    return 0;
}

void DFS(int CX[], long i)
{
    long j;
    CX[i] = 0;
    child[i].N = 0;

    for (j = 1; j <= N; j++)
        if (CX[j] && A_R[i][j])
        {
            child[i].N++;
            child[i].element[child[i].N] = j;
            DFS(CX, j);
        }
}

void DFS_delay(long i, List send[], long Fd[])
{
    long k, j, l, l2, l3, p, sl = 0;

    for (l = Fd[i] + 1; l <= MAXDELAY; l++)
    {
        if (sl == child[i].N)
            break;
        l2 = l % 20;
        l3 = l / 20 + 1;
        if (l2 == 0)
        {
            l2 = 20;
            l3 = l3 - 1;
        }

        if (n_color_slot[l2] != 0)
        {

            l3 = l3 % n_color_slot[l2];
            if (l3 == 0)
                l3 = n_color_slot[l2];

        /*for (k = 1; k <= send[i].N; k++)
        {
            j = send[i].element[k];
            if (j == l2 && (l3 % fre[i].element[k] == 0))
                break;
        }*/

            if ( inList(l2, send[i]) && (l3 == color_slot[l2][i]) )
                for (p = 1; p <= child[i].N; p++)
                {
                    j = child[i].element[p];
                    if (Fd[j] < 0)
                        if (inList(l2, grama[j]))
                        {
                            sl++;
                            Fd[j] = l;
                            if (child[j].N > 0)
                                DFS_delay(j, send, Fd);
                        }

                }
        }
    }

    if (l > MAXDELAY)
    {
        printf("exist error!");
        getch();
    }

}

void process()
{

    long i, j, k, l;
    double max, tmp;
    long truoc;
    double mauso;
    double r;
    double Fv[N_size + 5], Ft[N_size + 5];
    static Tree tT[N_size + 5];
    static Sche tB[N_size + 5];
    int CX[MAXN];
    int DD[MAXN];

    long delay = 0;
    long Fd[MAXN];

    start_time = clock();
    input();
    khoitao();

    max = -1;
    for (k = 1; k <= N_loop; k++)
    {
        mauso = 0;
        Ft[0] = 0;

        for (i = 1; i <= N_size; i++)
        {
            tmp = B[i].N * Es + T[i].M * Er;
            //Fv[i] = pow(maxE - tmp, alpha);
            //Fv[i] = 1.0 / pow(tmp, 1);
            Fv[i] = 1.0 / pow(tmp, alpha);
            if (max < Fv[i])
            {
                max = Fv[i];
                ganTree(T_KQ, T[i]);
                ganSche(B_KQ, B[i]);
            }
            mauso = mauso + Fv[i];
            Ft[i] = Fv[i] + Ft[i-1];
        }

        for (i = 1; i <= N_size; i++)
        {
            r = rand01();
            for (j = 1; j <= N_size; j++)
                if (Ft[j] + eps >= r * mauso)
                    break;

            ganTree(tT[i], T[j]);
            ganSche(tB[i], B[j]);
        }

        //quan the moi duoc lua chon la tT[i], tB[i]
        for (i = 1; i <= N_size; i++)   //gan lai cho T[i], B[i]
        {
            ganTree(T[i], tT[i]);
            ganSche(B[i], tB[i]);
        }

        truoc = 0;
        for (i = 1; i <= N_size; i++)
        {
            r = rand01();
            if (r <= Pc)
            {
                if (truoc == 0)
                    truoc = i;
                else    //lai ghep
                {
                    LGRand(tT[truoc], tB[truoc], tT[i], tB[i], T[truoc], B[truoc]);
                    LGWeight(tT[truoc], tB[truoc], tT[i], tB[i], T[i], B[i]);
                    truoc = 0;
                }
            }

        }   //end for i

        //dot bien
        for (i = 1; i <= N_size; i++)
        {
            r = rand01();
            if (r <= Pm)
                dotbien(T[i], B[i]);
        }
    }   //end for k

    for (i = 1; i <= N; i++)
        for (j = 1; j <= N; j++)
            A_R[i][j] = 0;

    for (i = 1; i <= T_KQ.M; i++)
    {
        //fprintf(fo, "%5ld%5ld\n", T_KQ.A[i], T_KQ.B[i]);
        A_R[T_KQ.A[i]][T_KQ.B[i]] = 1;
        A_R[T_KQ.B[i]][T_KQ.A[i]] = 1;
    }

    for (i = 1; i <= N; i++)
    {
        CX[i] = 1;
        send[i].N = 0;
        child[i].N = 0;
        DD[i] = 0;
    }
    for (i = 1; i <= T_KQ.M; i++)
    {
        DD[T_KQ.A[i]] = 1;
        DD[T_KQ.B[i]] = 1;
    }

    DFS(CX, S);

    for (i = 1; i <= N; i++)
    {
        int kt = 0;
        for (k = 1; k <= B_KQ.N; k++)
            if (B_KQ.Node[k] == i)
            {
                kt = 1;
                break;
            }

        if (kt)
        {
            //fprintf(fo, "%5ld:", i);
            for (k = 1; k <= B_KQ.N; k++)
            if (B_KQ.Node[k] == i)
            {
                //fprintf(fo, "%5ld", B_KQ.Slot[k]);
                send[i].N++;
                send[i].element[send[i].N] = B_KQ.Slot[k];
            }

            //fprintf(fo, "\n");
        }
    }

    //IN KQ
    fprintf(fo, "%ld\n", B_KQ.N);
    fprintf(fo, "%ld\n", B_KQ.N * Es + T_KQ.M * Er);


    //KHAI BAO CHO RIENG PHAN TINH DELAY
    long n_color, colored, i_color, vt, min, p;
    int need_color[MAXN];
    long color[MAXN], bac[MAXN];

    //NHO KHOI TAO send[j].N O TREN

    for (i = 1; i <= 20; i++)
    {

        n_color = 0;    //tong so luong dinh can phai to o slot i
        for (j = 1; j <= N; j++)
            need_color[j] = 0;

        for (j = 1; j <= N; j++)
            if (send[j].N > 0)     //la nut trong cua cay
            {
                if (inList(i, send[j]))
                {
                    need_color[j] = 1;
                    color[j] = 0;
                    n_color++;
                    for (p = 1; p <= N; p++)
                        X_conflict[j][p] = 0;   //ma tran ke cua do thi dung do
                }

            }

        for (j = 1; j <= N; j++)
            if (need_color[j])
                for (k = 1; k <= N; k++)
                    if (need_color[k])
                    {
                        if (A[j][k])                        //DIEU KIEN DUNG DO 1
                        {
                            X_conflict[j][k] = 1;
                            X_conflict[k][j] = 1;
                        }

                        for (l = 1; l <= adj[j].N; l++)     //DIEU KIEN DUNG DO 2
                        {
                            p = adj[j].element[l];
                            if (DD[p] && inList(i, grama[p]) && A[k][p])    //p thuoc cay ket qua va cung duoc ca j lan k truyen
                            {
                                X_conflict[j][k] = 1;
                                X_conflict[k][j] = 1;
                            }
                        }

                    }

        colored = 0;
        i_color = 0;
        while (colored < n_color)
        {
            i_color++;
            for (j = 1; j <= N; j++)
                if (need_color[j] && color[j] == 0)
                {
                    CX[j] = 1;
                }
                else
                    CX[j] = 0;
            while (1)   //chi to nhung dinh j co CX[j] = 1
            {
                for (j = 1; j <= N; j++)    //tinh bac
                    if (CX[j])
                        bac[j] = 0;

                for (j = 1; j <= N; j++)
                    if (CX[j])
                    {
                        for (k = 1; k <= N; k++)
                        {
                            if (CX[k] && X_conflict[j][k])
                                bac[j]++;
                        }
                    }
                min = MAXCOST;
                for (j = 1; j <= N; j++)
                    if (CX[j] && min > bac[j])
                    {
                        min = bac[j];
                        vt = j;
                    }
                if (min == MAXCOST)
                    break;

                color[vt] = i_color;
                colored++;
                CX[vt] = 0;
                for (k = 1; k <= N; k++)
                {
                    if (CX[k] && X_conflict[vt][k])
                        CX[k] = 0;
                }


            }   //end while (1)
        }   //end while

        n_color_slot[i] = i_color;
        for (j = 1; j <= N; j++)
            if (need_color[j])
            {
                color_slot[i][j] = color[j];
            }
    }

    for (i = 1; i <= N; i++)
        Fd[i] = -1;
    Fd[S] = 0;
    DFS_delay(S, send, Fd);
    for (i = 1; i <= N; i++)
        if (Fd[i] > delay)
            delay = Fd[i];

    fprintf(fo, "%5ld\n", delay);

    finish_time = clock();
    run_time = finish_time - start_time;
    run_time = run_time / CLOCKS_PER_SEC;
    fprintf(fo, "\n%20lf\n", run_time);

}

void ganList(List &A, List &B)
{
    long k;
    A.N = B.N;
    for (k = 1; k <= A.N; k++)
        A.element[k] = B.element[k];
}

void ganTree(Tree &T1, Tree &T2)
{
    long k;
    T1.M = T2.M;
    for (k = 1; k <= T1.M; k++)
    {
        T1.A[k] = T2.A[k];
        T1.B[k] = T2.B[k];
    }
}

void ganSche(Sche &A, Sche &B)
{
    long k;
    A.N = B.N;
    for (k = 1; k <= A.N; k++)
    {
        A.Node[k] = B.Node[k];
        A.Slot[k] = B.Slot[k];
    }
}

void dijRand(List adj[], List grama[], Tree &T3, Sche &B3)
//dijkstra voi trong so ngau nhien
{
    #define MAXWEIGHT 2000
    long i, j, k, l, min, sl, r;
    static long W[MAXN][MAXN];
    long D[MAXN];
    int CX[MAXN];
    long truoc[MAXN];

    for (i = 1; i <= N; i++)
        for (k = 1; k <= adj[i].N; k++)
        {
            j = adj[i].element[k];
            W[i][j] = rand() % MAXWEIGHT;
        }

    T3.M = 0;
    B3.N = 0;
    for (i = 1; i <= N; i++)
    {
        D[i] = MAXCOST;
        CX[i] = 1;
    }
    D[S] = 0;
    truoc[S] = 0;
    sl = K;
    while (sl > 0)
    {
        min = MAXCOST;
        for (j = 1; j <= N; j++)
            if (CX[j] && D[j] < min)
            {
                i = j;
                min = D[j];
            }

        CX[i] = 0;
        if (isTer[i])
            sl--;
        r = truoc[i];
        if (r > 0)
        {
            T3.M++;
            T3.A[T3.M] = r;
            T3.B[T3.M] = i;     //i la con cua r
            if (grama[i].N == 0)
            {
                printf("err %ld", i);
                getch();
            }
            l = grama[i].element[rand() % grama[i].N + 1];     //chon slot bat ki trong grama[i]
            for (k = 1; k <= B3.N; k++)     //kiem tra (r,l) da co trong lich truyen chua
                if (B3.Node[k] == r && B3.Slot[k] == l)
                    break;
            if (k > B3.N)
            {
                B3.N++;
                B3.Node[B3.N] = r;
                B3.Slot[B3.N] = l;
            }
        }

        for (k = 1; k <= adj[i].N; k++)
        {
            j = adj[i].element[k];
            if (CX[j] && D[j] > D[i] + W[i][j])
            {
                D[j] = D[i] + W[i][j];
                truoc[j] = i;
            }
        }

    }
}

void LGRand(Tree &T1, Sche &B1, Tree &T2, Sche &B2, Tree &T3, Sche &B3)
{
    long i, j, k, l, u, v, sl, r;
    long D[MAXN], truoc[MAXN];
    long min;
    static List adj2[MAXN], grama2[MAXN];

    for (i = 1; i <= N; i++)
    {
        adj2[i].N = 0;
        grama2[i].N = 0;
    }
    //hop canh
    for (i = 1; i <= T1.M; i++)
    {
        u = T1.A[i];
        v = T1.B[i];
        if (!inList(v, adj2[u]))
        {
            adj2[u].N++;
            adj2[u].element[adj2[u].N] = v;
        }
    }

    for (i = 1; i <= T2.M; i++)
    {
        u = T2.A[i];
        v = T2.B[i];
        if (!inList(v, adj2[u]))
        {
            adj2[u].N++;
            adj2[u].element[adj2[u].N] = v;
        }
    }

    //hop lich truyen
    //lich truyen cay 1
    for (i = 1; i <= N; i++)
    {
        for (k = 1; k <= T1.M; k++)
            if (T1.B[k] == i)
                break;

        if (k <= T1.M)      //dinh i la dinh toi cua cay T1
        {
            r = T1.A[k];    //r la cha cua i
            for (l = 1; l <= B1.N; l++)
                if (B1.Node[l] == r)
                {
                    j = B1.Slot[l];
                    if (inList(j, grama[i]) && !inList(j, grama2[i]))
                    {
                        grama2[i].N++;
                        grama2[i].element[grama2[i].N] = j;
                    }
                }
        }
    }   //end for i

    //lich truyen cay 2
    for (i = 1; i <= N; i++)
    {
        for (k = 1; k <= T2.M; k++)
            if (T2.B[k] == i)
                break;

        if (k <= T2.M)      //dinh i la dinh toi cua cay T2
        {
            r = T2.A[k];    //r la cha cua i
            for (l = 1; l <= B2.N; l++)
                if (B2.Node[l] == r)
                {
                    j = B2.Slot[l];
                    if (inList(j, grama[i]) && !inList(j, grama2[i]))
                    {
                        grama2[i].N++;
                        grama2[i].element[grama2[i].N] = j;
                    }
                }
        }
    }   //end for i
    dijRand(adj2, grama2, T3, B3);
}

void LGWeight(Tree &T1, Sche &B1, Tree &T2, Sche &B2, Tree &T3, Sche &B3)
{
    static List adj2[MAXN], grama2[MAXN], Child[MAXN];
    long i, j, k, l, r;
    long u, v;

    for (i = 1; i <= N; i++)
    {
        adj2[i].N = 0;
        grama2[i].N = 0;
    }
    //hop canh
    for (i = 1; i <= T1.M; i++)
    {
        u = T1.A[i];
        v = T1.B[i];
        if (!inList(v, adj2[u]))
        {
            adj2[u].N++;
            adj2[u].element[adj2[u].N] = v;
        }
    }

    for (i = 1; i <= T2.M; i++)
    {
        u = T2.A[i];
        v = T2.B[i];
        if (!inList(v, adj2[u]))
        {
            adj2[u].N++;
            adj2[u].element[adj2[u].N] = v;
        }
    }

    //hop lich truyen
    //lich truyen cay 1
    for (i = 1; i <= N; i++)
    {
        for (k = 1; k <= T1.M; k++)
            if (T1.B[k] == i)
                break;

        if (k <= T1.M)      //dinh i la dinh toi cua cay T1
        {
            r = T1.A[k];    //r la cha cua i
            for (l = 1; l <= B1.N; l++)
                if (B1.Node[l] == r)
                {
                    j = B1.Slot[l];
                    if (inList(j, grama[i]) && !inList(j, grama2[i]))
                    {
                        grama2[i].N++;
                        grama2[i].element[grama2[i].N] = j;
                    }
                }
        }
    }   //end for i

    //lich truyen cay 2
    for (i = 1; i <= N; i++)
    {
        for (k = 1; k <= T2.M; k++)
            if (T2.B[k] == i)
                break;

        if (k <= T2.M)      //dinh i la dinh toi cua cay T2
        {
            r = T2.A[k];    //r la cha cua i
            for (l = 1; l <= B2.N; l++)
                if (B2.Node[l] == r)
                {
                    j = B2.Slot[l];
                    if (inList(j, grama[i]) && !inList(j, grama2[i]))
                    {
                        grama2[i].N++;
                        grama2[i].element[grama2[i].N] = j;
                    }
                }
        }
    }   //end for i

    dij(adj2, grama2, ListTer, Child);
    T3.M = 0;
    B3.N = 0;
    for (i = 1; i <= N; i++)
        if (Child[i].N > 0)
        {
            for (k = 1; k <= Child[i].N; k++)
            {
                j = Child[i].element[k];
                T3.M++;
                T3.A[T3.M] = i;
                T3.B[T3.M] = j;
            }

            List temp;
            min_hit_set(grama2, Child[i], temp);

            for (k = 1; k <= temp.N; k++)
            {
                j = temp.element[k];
                B3.N++;
                B3.Node[B3.N] = i;
                B3.Slot[B3.N] = j;
            }
        }
}

void dotbien(Tree &Tc, Sche &Bc)
//CAN PHAN BIET adj (cua do thi input) vs adj2
{
    long i, j, k, l;
    long u, v, r;
    int CX[MAXN];   //CX[i] = 1 <=> chua thuoc cay Tc
    static List adj2[MAXN], grama2[MAXN], Child[MAXN];

    for (i = 1; i <= N; i++)
    {
        adj2[i].N = 0;
        grama2[i].N = 0;
        CX[i] = 1;
    }

    //lay canh
    for (i = 1; i <= Tc.M; i++)
    {
        u = Tc.A[i];
        v = Tc.B[i];
        CX[u] = 0;
        CX[v] = 0;
        if (!inList(v, adj2[u]))
        {
            adj2[u].N++;
            adj2[u].element[adj2[u].N] = v;
        }
    }

    //lay lich truyen
    for (i = 1; i <= N; i++)
        if (!CX[i])     //i thuoc Tc
        {
            for (k = 1; k <= Tc.M; k++)
            {
                j = Tc.B[k];
                if (j == i)
                    break;
            }

            if (k <= Tc.M)      //dinh i la dinh toi cua cay Tc
            {
                r = Tc.A[k];    //r la cha cua i
                for (l = 1; l <= Bc.N; l++)
                    if (Bc.Node[l] == r)
                    {
                        j = Bc.Slot[l];
                        if (inList(j, grama[i]) && !inList(j, grama2[i]))
                        {
                            grama2[i].N++;
                            grama2[i].element[grama2[i].N] = j;
                        }
                    }
            }


        }   //end for i

    //lay them cac dinh, canh ben ngoai
    for (i = 1; i <= N; i++)
        if (rand01() <= Ps)
        {
            ganList(adj2[i], adj[i]);
            ganList(grama2[i], grama[i]);
        }

    dij(adj2, grama2, ListTer, Child);
    Tc.M = 0;
    Bc.N = 0;
    for (i = 1; i <= N; i++)
        if (Child[i].N > 0)
        {
            for (k = 1; k <= Child[i].N; k++)
            {
                j = Child[i].element[k];
                Tc.M++;
                Tc.A[Tc.M] = i;
                Tc.B[Tc.M] = j;
            }

            List temp;
            min_hit_set(grama2, Child[i], temp);

            for (k = 1; k <= temp.N; k++)
            {
                j = temp.element[k];
                Bc.N++;
                Bc.Node[Bc.N] = i;
                Bc.Slot[Bc.N] = j;
            }
        }
}

void khoitao()
{
    long i;
    for (i = 1; i <= N_size; i++)
        dijRand(adj, grama, T[i], B[i]);
}

void input()
{
    long i, j, k;
    long u, v;

    start_time = clock();

    fi = fopen(FILEINP, "rt");
    fo = fopen(FILEOUT, "wt");
    fscanf(fi, "%ld%ld%ld%ld", &N, &M, &K, &S);  //so dinh, so canh, so terminal, dinh goc

    for (i = 1; i <= N; i++)
    {
        isTer[i] = 0;
        for (j = 1; j <= N; j++)
            A[i][j] = 0;
        adj[i].N = 0;
    }

    for (i = 1; i <= M; i++)
    {
        fscanf(fi, "%ld%ld", &u, &v);
        A[u][v] = 1;
        A[v][u] = 1;
        adj[u].N++;
        adj[v].N++;
        adj[u].element[adj[u].N] = v;
        adj[v].element[adj[v].N] = u;
    }

    for (i = 1; i <= N; i++)  //doc danh sach grama
    {
        fscanf(fi, "%ld", &grama[i].N);   //doc so luong phan tu
        for (j = 1; j <= grama[i].N; j++)
            fscanf(fi, "%ld", &grama[i].element[j]);
    }

    ListTer.N = 1;
    ListTer.element[1] = S;
    for (i = 1; i <= K; i++)    //doc so hieu cac Terminal
    {
        fscanf(fi, "%ld", &Ter[i]);
        isTer[Ter[i]] = 1;
        if (Ter[i] != S)
        {
            ListTer.N++;
            ListTer.element[ListTer.N] = Ter[i];
        }
    }
}

void dij(List adj[], List grama[], List &L, List Child[])
    //tap cac Terminal, nut nguon L.element[1]
    //KQ tra ra la Child cua cac nut
{
    long K, s;
    long i, j, k, l, min;
    long start;
    long D[MAXN], truoc[MAXN];
    int CX[MAXN], inTree[MAXN];
    List KQ, Tlist;
    static List adj2[MAXN];

    for (i = 1; i <= N; i++)
        adj2[i].N = 0;
    for (i = 1; i <= N; i++)
    {
        for (k = 1; k <= adj[i].N; k++)
        {
            j = adj[i].element[k];
            adj2[j].N++;
            adj2[j].element[adj2[j].N] = i;
        }
    }
    for (i = 1; i <= N; i++)
        ganList(adj[i], adj2[i]);

    K = L.N;
    s = L.element[1];

    for (i = 1; i <= N; i++)
    {
        inTree[i] = 0;
        Child[i].N = 0;
    }
    inTree[s] = 1;

    for (k = 2; k <= K; k++)
    {
        start = L.element[k];
        for (i = 1; i <= N; i++)
        {
            D[i] = MAXCOST;
            CX[i] = 1;
        }
        D[start] = 0;
        truoc[start] = 0;

        while (1)
        {
                min = MAXCOST;
                for (j = 1; j <= N; j++)
                    if (CX[j] && min > D[j])
                    {
                        min = D[j];
                        i = j;
                    }
                CX[i] = 0;
                if (min == MAXCOST)
                {
                    printf("sai roi");
                    getch();
                }
                if (inTree[i])
                    break;

                for (l = 1; l <= adj[i].N; l++)
                {
                    j = adj[i].element[l];
                    if (CX[j])
                    {
                        if (inTree[j])
                        {
                            ganList(Tlist, Child[j]);
                            Tlist.N++;
                            Tlist.element[Tlist.N] = i;

                            long tmp1 = min_hit_set(grama, Tlist, KQ);
                            long tmp2 = min_hit_set(grama, Child[j], KQ);
                            if (D[j] > D[i] + (tmp1 - tmp2) * Es + Er)
                            {
                                D[j] = D[i] + (tmp1 - tmp2) * Es + Er;
                                truoc[j] = i;
                            }

                        }
                        else
                        {
                            if (D[j] > D[i] + Es + Er)
                            {
                                D[j] = D[i] + Es + Er;
                                truoc[j] = i;
                            }
                        }

                    }
                }   //end for l

        }   //end while (1)
        while (truoc[i] != 0)
        {
            j = truoc[i];
            Child[i].N = Child[i].N + 1;
            Child[i].element[Child[i].N] = j;
            i = truoc[i];
            inTree[i] = 1;
        }
    }   //end for k
}

int inList(long x, List &L)
{
    long i;
    for (i = 1; i <= L.N; i++)
        if (L.element[i] == x)
            break;
    if (i <= L.N)
        return 1;
    else
        return 0;
}

long min_hit_set(List grama[], List &Child, List &KQ)  //tra ra KQ la cac SLOT
{
    long i, j, k, l, vt;
    long max, tong, max_tong;
    long N_hit;
    int CX[MAXSLOT];
    int hit[MAXN];

    max = 0;
    for (i = 1; i <= Child.N; i++)
    {
        j = Child.element[i];
        for (k = 1; k <= grama[j].N; k++)
            if (max < grama[j].element[k])
                max = grama[j].element[k];
    }

    for (i = 1; i <= max; i++)
        CX[i] = 0;

    for (i = 1; i <= Child.N; i++)
    {
        j = Child.element[i];
        for (k = 1; k <= grama[j].N; k++)
            CX[grama[j].element[k]] = 1;
    }

    N_hit = 0;
    KQ.N = 0;
    for (i = 1; i <= Child.N; i++)
        hit[i] = 0;

    while (N_hit < Child.N)
    {
        max_tong = 0;
        for (i = 1; i <= max; i++)
            if (CX[i])
            {
                tong = 0;
                for (k = 1; k <= Child.N; k++)
                {
                    j = Child.element[k];
                    if (hit[k] == 0 && inList(i, grama[j]))
                    {
                        tong++;
                    }
                }

                if (max_tong < tong)
                {
                    max_tong = tong;
                    vt = i;
                }
            }

        KQ.N++;
        KQ.element[KQ.N] = vt;
        CX[vt] = 0;

        for (k = 1; k <= Child.N; k++)
        {
            j = Child.element[k];
            if (hit[k] == 0 && inList(vt, grama[j]))
            {
                N_hit++;
                hit[k] = 1;
            }
        }
    }       //end while
    return KQ.N;
}

double rand01()
{
    return 1.0 * rand() / RAND_MAX;
}
