#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <conio.h>

#define MAXL 2000
#define MAXADJ 2000
#define MAXN 2000
#define MAXSLOT 35
#define MAXCOST 10000000
#define NLOOP 1
#define MAXDELAY 10000000

FILE *fi, *fo;

char *FILEINP, *FILEOUT;
struct List
{
    long N;
    long element[MAXL];
};

//khai bao
long Es = 100, Er = 15;
clock_t start_time, finish_time;
double run_time;

long N, M, K, S, M_tree, M_KQ;
List adj[MAXN], grama[MAXN];
long A_tree[MAXN], B_tree[MAXN], A_KQ[MAXN], B_KQ[MAXN];

int A[MAXN][MAXN], A_tmp[MAXN][MAXN];
long Ter[MAXN];
List Child[MAXN], Child_KQ[MAXN], Child_LS[MAXN], Ter_LS, L_tmp;
List send[MAXN];
int inTree[MAXN];
long source[MAXN];

int X_conflict[MAXN][MAXN];
long n_color_slot[MAXSLOT], color_slot[MAXSLOT][MAXN];

//nguyen mau ham

void input();
void find_path(long []);
long min_hit_set(List, List &);
void hoanvi(long &, long &);
void dij(List, List[]);
void DFS(int [], List [], long, long );
void DFS(int [], long );
long total_energy(long, List[]);
int inList(long , List );

void ganList(List &A, List &B)
{
    long k;
    A.N = B.N;
    for (k = 1; k <= A.N; k++)
        A.element[k] = B.element[k];
}

void DFS_delay(long i, List send[], long Fd[])
{
    long k, j, l, l2, l3, p, sl = 0;

    for (l = Fd[i] + 1; l <= MAXDELAY; l++)
    {
        if (sl == Child_KQ[i].N)
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
                for (p = 1; p <= Child_KQ[i].N; p++)
                {
                    j = Child_KQ[i].element[p];
                    if (Fd[j] < 0)
                        if (inList(l2, grama[j]))
                        {
                            sl++;
                            Fd[j] = l;
                            if (Child_KQ[j].N > 0)
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
    long sl = 0;
    long i, j, k, l, min, total, loop;
    long tmp1, tmp2;
    long trans;
    int CX[MAXN];
    int DD[MAXN];

    long delay = 0;
    long Fd[MAXN];

    input();
    for (i = 1; i <= K; i++)
        if (Ter[i] != S)
        {
            sl++;
            source[sl] = Ter[i];
        }

    min = MAXCOST;  //kq total energy
    for (loop = 1; loop <= NLOOP; loop++)
    {
        //for (j = 1; j < K; j++)
            //hoanvi(source[j], source[rand() % (K - 1) + 1]);

        find_path(source);

        //them LOCAL SEARCH
        for (i = 1; i <= N; i++)   //duyet cac con cua nut S
        if (inTree[i])  //i thuoc cay hien tai
        {
            for (j = 1; j <= N; j++)
                for (k = 1; k <= N; k++)
                    A_tmp[j][k] = 0;    //ma tran ke cua cay moi

            //i = Child[S].element[k];
            //A_tmp[i][S] = 1;
            //A_tmp[S][i] = 1;
            tmp1 = total_energy(i, Child);

            for (j = 1; j <= N; j++)
                CX[j] = 1;
            DFS(CX, Child, S, i);  //dua ra do thi = do thi goc - cay con goc i

            Ter_LS.N = 1;   //tap Terminal cua cay con goc i
            Ter_LS.element[1] = i;
            for (l = 1; l <= K; l++)
                if (CX[Ter[l]] && Ter[l] != i)
                {
                    Ter_LS.N++;
                    Ter_LS.element[Ter_LS.N] = Ter[l];
                }

            if (Ter_LS.N > 2)   //goc i co nhieu hon 2 Terminal
            {
                //hoan vi ngau nhien tap Ter
                for (j = 2; j <= Ter_LS.N; j++)
                    hoanvi(Ter_LS.element[j], Ter_LS.element[rand() % (Ter_LS.N - 1) + 2]);
                dij(Ter_LS, Child_LS);

                tmp2 = total_energy(i, Child_LS);
                if (tmp1 > tmp2)
                {
                    for (j = 1; j <= N; j++)
                        CX[j] = 1;
                    DFS(CX, Child_LS, i, -1);  //lay ra do thi goc i moi

                    M_tree = 0;
                    for (j = 1; j <= N; j++)
                    {
                        Child[j].N = 0;
                        inTree[j] = 0;
                    }
                    for (j = 1; j <= N; j++)
                        CX[j] = 1;
                    DFS(CX, S); //dua cay ket qua vao Child, inTree, M_tree
                }

            }

        }
        //ket thuc LOCAL SEARCH

        total = 0;
        for (j = 1; j <= N; j++)
            if (inTree[j] && Child[j].N > 0)
            {
                List tmp;
                min_hit_set(Child[j], tmp);
                total = total + tmp.N * Es + Child[j].N * Er;
            }

        if (min > total)
        {
            min = total;
            M_KQ = M_tree;
            for (j = 1; j <= N; j++)
            {
                Child_KQ[j].N = 0;
                DD[j] = 0;

                if (inTree[j] && Child[j].N > 0)
                {
                    Child_KQ[j] = Child[j];
                }
            }

            for (j = 1; j <= M_KQ; j++)
            {
                A_KQ[j] = A_tree[j];
                B_KQ[j] = B_tree[j];
                DD[A_KQ[j]] = 1;
                DD[B_KQ[j]] = 1;
            }
        }

    }   //end for loop

    //___________IN KQ____________
    for (i = 1; i <= M_KQ; i++)
    {
        //fprintf(fo, "%5ld %5ld\n", A_KQ[i], B_KQ[i]);
    }
    trans = 0;
    for (i = 1; i <= N; i++)
        if (Child_KQ[i].N > 0)
        {
            List tmp;
            min_hit_set(Child_KQ[i], tmp);
            trans = trans + tmp.N;
            ganList(send[i], tmp);

            /*fprintf(fo, "%5ld:", i);
            for (j = 1; j <= tmp.N; j++)
                fprintf(fo, "%5ld", tmp.element[j]);
            fprintf(fo, "\n");*/
        }
        else
            send[i].N = 0;

    fprintf(fo, "%ld\n", trans);
    fprintf(fo, "%ld\n", min);

    //KHAI BAO CHO RIENG PHAN TINH DELAY
    long n_color, colored, i_color, vt, p;
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

//int main()
int main(int Nts, char **arg)
{

    srand(time(NULL));
    FILEINP = arg[1];
    FILEOUT = arg[2];
    //FILEINP = "WSN03.inp";
    //FILEOUT = "WSN03.outD";
    process();
    fclose(fi);
    fclose(fo);
    return 0;
}

void hoanvi(long &x, long &y)
{
    long tmp;
    tmp = x;
    x = y;
    y = tmp;
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

    for (i = 1; i <= K; i++)    //doc so hieu cac Terminal
    {
        fscanf(fi, "%ld", &Ter[i]);
    }
}

int inList(long x, List L)
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

long min_hit_set(List Child, List &KQ)  //tra ra KQ la cac SLOT
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

void find_path(long source[])   //tim duong di tu cac terminal Source den S
{
    long i, j, k, t, start;
    List Tlist, KQ;
    long D[MAXN];
    int CX[MAXN];
    long truoc[MAXN];
    long min;

    M_tree = 0;
    for (i = 1; i <= N; i++)
    {
        inTree[i] = 0;
        Child[i].N = 0;
    }

    inTree[S] = 1;


    for (t = 1; t < K; t++)     //xet cac terminal
    {
        start = source[t];
        for (i = 1; i <= N; i++)
        {
            D[i] = MAXCOST;
            CX[i] = 1;
        }

        truoc[start] = 0;
        D[start] = 0;

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
            if (inTree[i])  //i la dinh thuoc cay truoc
                break;

            for (k = 1; k <= adj[i].N; k++)
            {
                j = adj[i].element[k];
                if (CX[j])
                {
                    if (inTree[j])
                    {
                        Tlist = Child[j];
                        Tlist.N = Child[j].N + 1;
                        Tlist.element[Tlist.N] = i;

                        long tmp1 = min_hit_set(Tlist, KQ);
                        long tmp2 = min_hit_set(Child[j], KQ);
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

                }   //end if

            }   //end for j
        }       //end while


        //__________cap nhat cay_________

        while (truoc[i] != 0)
        {
            j = truoc[i];
            Child[i].N++;
            Child[i].element[Child[i].N] = j;
            M_tree++;
            A_tree[M_tree] = i;
            B_tree[M_tree] = j;
            i = truoc[i];
            inTree[i] = 1;
        }

    }   //end for t



}

void dij(List L, List Child[])    //tap cac Terminal, nut nguon L.element[1]
                                //KQ tra ra la Child cua cac nut
{
    long K, s;
    long i, j, k, l, min;
    long start;
    long D[MAXN], truoc[MAXN];
    int CX[MAXN], inTree[MAXN];
    List KQ, Tlist;

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

                if (inTree[i])
                    break;

                for (l = 1; l <= adj[i].N; l++)
                {
                    j = adj[i].element[l];
                    if (CX[j])
                    {
                        if (inTree[j])
                        {
                            Tlist = Child[j];
                            Tlist.N++;
                            Tlist.element[Tlist.N] = i;

                            long tmp1 = min_hit_set(Tlist, KQ);
                            long tmp2 = min_hit_set(Child[j], KQ);
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

void DFS(int CX[], List Child[], long i, long r)    //DFS tranh cay con goc r
{
    long k, l, j;

    CX[i] = 0;
    for (k = 1; k <= Child[i].N; k++)
    {
        j = Child[i].element[k];
        if (CX[j])
        {
            A_tmp[i][j] = 1;
            A_tmp[j][i] = 1;
            if (j != r)
                DFS(CX, Child, j, r);
        }
    }

}

void DFS(int CX[], long i)  //dua cay ket qua vao Child, inTree, M_tree
{
    long k, l, j;

    CX[i] = 0;
    inTree[i] = 1;
    for (j = 1; j <= N; j++)
        if (CX[j] && A_tmp[i][j])
        {
            M_tree++;
            A_tree[M_tree] = i;
            B_tree[M_tree] = j;
            Child[i].N++;
            Child[i].element[Child[i].N] = j;
            DFS(CX, j);
        }

}

long total_energy(long i, List Child[])
{
    long t, j, k;
    if (Child[i].N > 0)
    {
        t = 0;
        for (k = 1; k <= Child[i].N; k++)
        {
            j = Child[i].element[k];
            t = t + total_energy(j, Child);
        }
        return t + min_hit_set(Child[i], L_tmp) * Es + Child[i].N * Er;
    }
    else return 0;
}
