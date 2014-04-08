#include <stdio.h>
#include <time.h>
#include <conio.h>

#define MAXL 6305
#define MAXADJ 6305
#define MAXN 6305
#define MAXSLOT 35
#define MAXCOST 10000000
#define MAXNE 6305
#define MAXDELAY 10000000

FILE *fi, *fo;
//char *FILEINP ="WSN42.inp", *FILEOUT = "WSN42.outH";
char *FILEINP ="", *FILEOUT = "";

clock_t start_time, finish_time;
double run_time;
struct List
{
    long N;
    long element[MAXL];
};

struct Node
{
    int type;
    //List grama;
};

//khai bao
long Es = 100, Er = 15;
long N, M, K, S;
long N_E, M_E, K_E;
List adj[MAXN], grama[MAXN];
List adj_E[MAXNE];
List send[MAXN];//, fre[MAXN];

int X_conflict[MAXN][MAXN];
int A[MAXN][MAXN];
int A_E[MAXNE][MAXNE];
int A_R[MAXN][MAXN];    //canh cua cay result
//long N_child[MAXN], child[MAXN];
List child[MAXN];

long A_S[MAXNE][MAXNE];
//long P[MAXN][MAXN];

long Ter[MAXN];
long C[MAXN][MAXSLOT];
long nuclear[MAXN*MAXSLOT + MAXN], slot[MAXN*MAXSLOT + MAXN];
long sl;    //so luong satelli duoc chon
long Satelli[MAXN];
int X[MAXNE][MAXNE];  //chi dung de loai bo dinh thua
int use[MAXNE];     //danh dau cac satellite co mat trong cay steiner cuoi
long M_S;   //so luong canh cua cay steiner satellite
long A1[MAXNE], A2[MAXNE];    //luu tru canh cua cay steiner satellite

long n_color_slot[MAXSLOT], color_slot[MAXSLOT][MAXN];

//nguyen mau ham
void input();
void extend_graph();
void select_satellite();
void steiner_tree();
void result_tree();

void process()
{
    input();
    extend_graph();
    select_satellite();
    steiner_tree();
    result_tree();
}

int main(int Nts, char **arg)
//int main()
{
    FILEINP = arg[1];
    FILEOUT = arg[2];
    process();
    fclose(fi);
    fclose(fo);
    return 0;
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

void extend_graph()
{
    //int CX[MAXN][MAXSLOT];
    long i, j, k, l, t, a, b, u, v;

    N_E = N;
    //M_E = 0;
    for (i = 1; i <= N; i++)
    {
        nuclear[i] = i;
        for (j = 1; j < MAXSLOT; j++)
            C[i][j] = 0;   //chua co dinh (i,j)

        for (k = 1; k <= adj[i].N; k++)
        {
            t = adj[i].element[k];  //t la mot dinh ke voi i
            for (l = 1; l <= grama[t].N; l++)
            {
                j = grama[t].element[l];
                if (C[i][j] == 0)   //them dinh (i, j) vao do thi MR
                {
                    N_E++;
                    nuclear[N_E] = i;
                    slot[N_E] = j;
                    C[i][j] = N_E;
                }
            }
        }

    }   //end for i

    //xay dung canh cho do thi mo rong
    for (i = 1; i <= N_E; i++)
        for (j = 1; j <= N_E; j++)
            A_E[i][j] = 0;

    //xay dung do thi con day du
    for (i = 1; i <= N_E; i++)
    {
        adj_E[i].N = 0;
    }

    for (i = 1; i <= N_E; i++)
        for (j = 1; j <= N_E; j++)
            if (i != j && nuclear[i] == nuclear[j])
            {
                A_E[i][j] = 1;
                A_E[j][i] = 1;
                adj_E[i].N++;
                adj_E[j].N++;
                adj_E[i].element[adj_E[i].N] = j;
                adj_E[j].element[adj_E[j].N] = i;
            }

    for (u = 1; u <= N; u++)
        for (v = 1; v <= N; v++)
            if (A[u][v])
            {
                //noi (v,i) voi u
                for (k = 1; k <= grama[u].N; k++)
                {
                    i = grama[u].element[k];
                    b = C[v][i];
                    A_E[u][b] = 1;
                    A_E[b][u] = 1;
                    adj_E[u].N++;
                    adj_E[b].N++;
                    adj_E[u].element[adj_E[u].N] = b;
                    adj_E[b].element[adj_E[b].N] = u;
                }

                //noi (u,j) voi v
                for (l = 1; l <= grama[v].N; l++)
                {
                    j = grama[v].element[l];
                    a = C[u][j];
                    A_E[a][v] = 1;
                    A_E[v][a] = 1;
                    adj_E[a].N++;
                    adj_E[v].N++;
                    adj_E[a].element[adj_E[a].N] = v;
                    adj_E[v].element[adj_E[v].N] = a;
                }

                //noi (u,j) voi (v,i)
                for (k = 1; k <= grama[u].N; k++)
                    for (l = 1; l <= grama[v].N; l++)
                    {
                        i = grama[u].element[k];
                        j = grama[v].element[l];
                        a = C[u][j];
                        b = C[v][i];

                        //noi canh
                        A_E[a][b] = 1;
                        A_E[b][a] = 1;
                        adj_E[a].N++;
                        adj_E[b].N++;
                        adj_E[a].element[adj_E[a].N] = b;
                        adj_E[b].element[adj_E[b].N] = a;
                    }

            }
}

void select_satellite()
{
    long i, j, k, l;
    long max, tong, vt;
    int CX[MAXN];
    long N_hit;

    N_hit = 0;
    sl = 0;
    for (i = 1; i <= K; i++)
        CX[i] = 1;      //chua cham vao terminal thu i

    while (N_hit < K)
    {
        max = 0;
        for (i = N + 1; i <= N_E; i++)  //lua chon cac ve tinh
        {
            tong = 0;
            for (k = 1; k <= K; k++)
            {
                j = Ter[k];     //ve tinh thu k
                if (A_E[i][j] && CX[k])     //chua xet ve tinh thu k
                {
                    tong++;
                }
            }
            if (max < tong)
            {
                vt = i;
                max = tong;
            }
        }   //end for i

        sl++;
        Satelli[sl] = vt;

        for (k = 1; k <= K; k++)    //danh dau cac terminal lien ke voi Satelli duoc chon
        {
                j = Ter[k];
                if (A_E[vt][j] && CX[k])
                {
                    N_hit++;
                    CX[k] = 0;
                }
        }

    }   //end while
}

void find_path(long j, long truoc[])
{
    long i;
    while (truoc[j] > 0)
    {
        i = truoc[j];
        A_S[i][j] = 1;
        A_S[j][i] = 1;
        j = truoc[j];
    }
}

void path_tree(int CX[], long s, long t)
{
    //int CX[MAXNE];
    long d, c, i, j, k;
    long Q[MAXNE];

    d = 1;
    c = 1;
    Q[1] = s;
    use[s] = 1;

    for (j = N + 1; j <= N_E; j++)
        CX[j] = 1;
    CX[s] = 0;

    while (d <= c)
    {
        i = Q[d];
        for (j = N + 1; j <= N_E; j++)
            if (CX[j] && X[i][j] > 0)
            {
                c++;
                Q[c] = j;
                CX[j] = 0;
                use[j] = 1;
                if (j == t)
                    return;
            }
        d++;
    }
}

void dij(long s, long t, long F[], long truoc[])
//tim duong di tu s den t tren do thi mo rong, kq tra ve trong mang F[], truoc[]
//neu t = 0 thi la tim duong di tu s den tat ca cac nut Satelli
{
    int CX[MAXNE], is_Satelli[MAXNE];
    int lap;
    long i, j, k, dau, cuoi, tong;
    long Q[MAXNE];

    for (i = N + 1; i <= N_E; i++)
    {
        F[i] = MAXCOST;
        is_Satelli[i] = 0;
        CX[i] = 1;
    }

    for (i = 1; i <= sl; i++)
        is_Satelli[Satelli[i]] = 1;

    dau = 1;    cuoi = 1;
    Q[1] = s;   CX[s] = 0;  truoc[s] = 0;
    F[s] = 0;   lap = 1;
    if (is_Satelli[s])
        tong = sl - 1;
    else
        tong = sl;
    while (lap && dau <= cuoi)
    {
        i = Q[dau];
        for (k = 1; k <= adj_E[i].N; k++)
        {
            j = adj_E[i].element[k];
            if (CX[j] && j >= N + 1)
            {
                cuoi++;
                Q[cuoi] = j;
                CX[j] = 0;
                F[j] = F[i] + 1;
                truoc[j] = i;
                if (t == 0)
                {
                    if (is_Satelli[j])
                    {
                        tong--;
                        if (tong == 0)
                        {
                            lap = 0;
                            break;
                        }
                    }
                }
                else
                    if (j == t)
                    {
                        lap = 0;
                        break;
                    }
            }
        }
        dau++;
    }

}
void steiner_tree()
{
    long i, j, k, l;
    long F[MAXNE], truoc[MAXNE];
    static long W[MAXN][MAXN];
    long prev[MAXNE], D[MAXNE], min;
    int CX[MAXNE];

    for (k = 1; k <= sl; k++)
    {
        i = Satelli[k];
        dij(i, 0, F, truoc);
        for (l = 1; l <= sl; l++)
            if (l != k)
            {
                j = Satelli[l];
                W[k][l] = F[j];
                W[l][k] = W[k][l];

            }
            else
            {
                W[k][l] = MAXCOST;
                W[l][k] = W[k][l];
            }
    }

    prev[1] = 0;
    for (i = 1; i <= sl; i++)
    {
        D[i] = MAXCOST;
        CX[i] = 1;
    }

    D[1] = 0;

    for (k = 1; k <= sl; k++)
    {
        min = MAXCOST;
        for (j = 1; j <= sl; j++)
            if (CX[j] && D[j] < min)
            {
                min = D[j];
                i = j;
            }

        CX[i] = 0;
        for (j = 1; j <= sl; j++)
            if (CX[j] && D[j] > W[i][j])
            {
                D[j] = W[i][j];
                prev[j] = i;
            }
    }

    for (i = N + 1; i <= N_E; i++)
        for (j = N + 1; j <= N_E; j++)
            A_S[i][j] = MAXCOST;
    //liet ke cac duong di
    for (k = 2; k <= sl; k++)
    {
        i = Satelli[k];
        j = Satelli[prev[k]];
        dij(i, j, F, truoc);
        //find_path(Satelli[k], Satelli[prev[k]]);
        find_path(j, truoc);
    }

    //tim cay khung lan 2, A_S[i][j] la trong so
    for (i = N + 1; i <= N_E; i++)
    {
        D[i] = MAXCOST;
        CX[i] = 1;
        prev[i] = 0;
    }
    D[Satelli[1]] = 0;
    prev[Satelli[1]] = 0;

    while (1)
    {
        min = MAXCOST;
        for (j = N + 1; j <= N_E; j++)
            if (CX[j] && D[j] < min)
            {
                min = D[j];
                i = j;
            }
        if (min == MAXCOST)
            break;
        CX[i] = 0;
        for (j = N + 1; j <= N_E; j++)
            if (CX[j] && D[j] > A_S[i][j])
            {
                D[j] = A_S[i][j];
                prev[j] = i;
            }

    }
    // loai bo cac dinh thua
    for (i = N + 1; i <= N_E; i++)
        for (j = N + 1; j <= N_E; j++)
            X[i][j] = 0;    //X bieu thi canh cua cay thu duoc

    for (i = N + 1; i <= N_E; i++)
    {
        if (prev[i] > 0)
        {
            X[i][prev[i]] = 1;
            X[prev[i]][i] = 1;
        }
    }

    //loai bo bang cach tim duong di tu satelli[1] den cac satelli khac
    for (i = 0; i <= N_E; i++)
        use[i] =  0;
    for (i = 2; i <= sl; i++)
    {
        for (j = N + 1; j <= N_E; j++)
            CX[j] = 1;
        path_tree(CX, Satelli[1], Satelli[i]);

    }
    M_S = 0;    //so luong canh cua cay steiner cuoi cung (kq)
    for (i = N + 1; i <= N_E; i++)
    {
        j = prev[i];
        if (j > 0 && use[i] && use[j])
        {
            M_S++;
            A1[M_S] = i;
            A2[M_S] = j;
        }
    }
}

void DFS(int DD[], int CX[], long i)
{
    long j;
    CX[i] = 0;
    child[i].N = 0;

    for (j = 1; j <= N; j++)
        if (DD[j] && CX[j] && A_R[i][j])
        {
            //fprintf(fo, "%5ld %5ld\n", i, j);
            child[i].N++;
            child[i].element[child[i].N] = j;
            DFS(DD, CX, j);
        }
}

void ganList(List &A, List &B)
{
    long k;
    A.N = B.N;
    for (k = 1; k <= A.N; k++)
        A.element[k] = B.element[k];
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
void result_tree()
{
    int DD[MAXN];   //cac dinh tren do thi goc
    int CX[MAXN];
    long i, j, k, l, p;
    long total = 0;
    long totalE = 0;
    long slnhan = 0;
    long delay = 0;
    long Fd[MAXN];

    for (i = 1; i <= N; i++)
    {
        for (j = 1; j <= N; j++)
        {
            A_R[i][j] = 0;  //A_R bieu dien canh cua cay ket qua
        }
        DD[i] = 0;
    }

    for (k = 1; k <= M_S; k++)
    {
        i = A1[k];
        j = A2[k];
        A_R[nuclear[i]][nuclear[j]] = 1;
        A_R[nuclear[j]][nuclear[i]] = 1;
        DD[nuclear[i]] = 1;     //DD - danh dau nhung dinh thuoc cay cuoi cung
        DD[nuclear[j]] = 1;
    }

    for (k = 1; k <= K; k++)    //noi cac Terminal den dinh lien ke
    {
        i = Ter[k];
        if (DD[i] == 0)
        {
            for (j = N + 1; j <= N_E; j++)
                if (use[j] && A_E[i][j])
                {
                    A_R[i][nuclear[j]] = 1;
                    A_R[nuclear[j]][i] = 1;
                    break;  //noi 1 canh duy nhat den nut terminal i
                }
            DD[i] = 1;
        }
    }

    for (i = 1; i <= N; i++)
    {
        CX[i] = 1;
        child[i].N = 0;
        send[i].N = 0;
    }

    DFS(DD, CX, S);     //chi xet den dinh i co DD[i] = 1
    total = 0;  //so luong slot truyen
    totalE = 0; //tong nang luong truyen

    for (i = 1; i <= N; i++)
        if (DD[i] && child[i].N > 0)    //DD[i] = 1 <=> i thuoc cay cuoi cung
        {
            List temp;
            if (child[i].N == 1)    //i thuoc d1(T)
            {
                j = child[i].element[1];
                temp.N = 1;
                temp.element[1] = grama[j].element[1];
            }
            else    //i thuoc d+(T)
            {
                temp.N = 0;
                for (j = N + 1; j <= N_E; j++)
                    if (use[j] && nuclear[j] == i)  //dua ra F(i)
                    {
                        temp.N++;
                        temp.element[temp.N] = slot[j];
                    }
            }

            total = total + temp.N;
            totalE = totalE + temp.N * Es;
            ganList(send[i], temp);

            /*fprintf(fo, "%5ld:", i);
            for (j = 1; j <= temp.N; j++)
                fprintf(fo, "%5ld", temp.element[j]);
            fprintf(fo, "\n");*/
        }

    slnhan = 0;
    for (i = 1; i <= N; i++)
        if (DD[i])
        {
            slnhan = slnhan + 1;
        }


    slnhan = slnhan - 1;    //ngoai tru dinh S
    totalE = totalE + slnhan * Er;
    fprintf(fo, "%5ld\n", total);
    fprintf(fo, "%5ld\n", totalE);

    /*for (i = 1; i <= N; i++)
        Fd[i] = -1;
    Fd[S] = 0;
    DFS_delay(S, send, Fd);
    for (i = 1; i <= N; i++)
        if (Fd[i] > delay)
            delay = Fd[i];*/

    //KHAI BAO CHO RIENG PHAN TINH DELAY
    long n_color, colored, i_color, min, vt;
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
    fprintf(fo, "%20lf\n", run_time);
}
