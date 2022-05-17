#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include<bits/stdc++.h>
#include <omp.h>


using namespace std;
using std::max;


// ----- Pesos definidos ----- //
#define WMAT 2
#define WMIS -1
#define WGAP -1


// ----- Structs ----- //
struct item {
    int item_score;
    vector<string> seq_a;
    vector<string> seq_b;
};

struct item_sw{
    int valor;
    int salto_i;
    int salto_j;
};



// ----- Declarações de funções ----- //
int m,n;
int w(string a, string b);
void gera_subseq(vector<string> seq, int start_point, int end_point, vector<vector<string>>& matriz_subseq);
item smith_waterman(vector<string> a, vector<string> b);



int main() {
    double init_time, final_time;
    init_time = omp_get_wtime();


    int n=0;
    vector<string> a;
    vector<string> b;
    
 // Lendo tamanho das sequências   
    cin >> n >> m;

// Lendo sequências iniciais
    string base;
    cin >> base;
    for(int i = 0; i < n; i++){
        a.push_back({base[i]});
    }
    cin >> base;
    for(int i = 0; i < m; i++){
        b.push_back({base[i]});
    }

    vector<vector<string>> subseqs_a;
    vector<vector<string>> subseqs_b;
    gera_subseq(a, 0, 1, subseqs_a);
    gera_subseq(b, 0, 1, subseqs_b);

// ----- Obtendo o score e as sequências alinhadas ----- //
    
    item melhor, sw_atual;    
    vector<item> combinations;   
    int i=0;
    for (auto& sub_a : subseqs_a){
        for (auto& sub_b : subseqs_b){
            combinations.push_back({i,sub_a, sub_b});
            i+=1;
        }
    }

    vector<item> resultados((int)combinations.size());
    #pragma omp parallel for shared(resultados) 
    for (auto& el : combinations){
        resultados[el.item_score] = smith_waterman(el.seq_a, el.seq_b);
    }

    int melhor_valor=-1;
    for (int i=0; i<(int)resultados.size(); i++){
        if (resultados[i].item_score > melhor_valor){
            melhor_valor = resultados[i].item_score;
            melhor = resultados[i];
        }
     }

    

// ----- Imprimindo Score obtido ----- //
    cout << "----- Score -----" << endl;
    cout << "max_score: "<< melhor.item_score << endl << endl;

    cout << " ----- Melhor subsequência de A -----" << endl;
    for (auto& el : melhor.seq_a ){
        cout << el << " ";
    }
    cout << endl;

    cout << " ----- Subsequência B correspondente -----" << endl;
    for (auto& el : melhor.seq_b){
        cout << el << " ";
    }
    cout << endl;



    final_time = omp_get_wtime() - init_time;
    cout << "tempo: " << final_time << endl;

    // cout << melhor.item_score ;

     return 0;

}


// ----- Função que gera todas as subsequências possíveis ----- //
void gera_subseq(vector<string> seq, int start_point, int end_point, vector<vector<string>>& matriz_subseq){
    if (end_point > (int)seq.size())
      return;
    else if (start_point > end_point){
        gera_subseq(seq, 0, end_point+1, matriz_subseq);
    }else{
        if (start_point != end_point){
            vector<string> subseq;
            for (int j=start_point; j<end_point; j++){
                subseq.push_back(seq[j]);
            }
            matriz_subseq.push_back(subseq);
        }
        gera_subseq(seq, start_point+1, end_point, matriz_subseq);
    }
}



// ----- Função que retorna os pesos do cálculo do score de Smith Waterman ----- //
int w(string a, string b){
    if (a == b)
        return WMAT;  //match
    else if (a != b)
        return WMIS;  // mismatch
   else
        return WGAP;  // gap

}


// ----- Função que acha o máximo dentre as opções de diagonal, inserção ou deleção para o cálculo de Smith Waterman ----- //
item_sw find_max(int diagonal, int delecao, int insercao){
    item_sw maior;
    int max_score_local=0;
    max_score_local = max({0, diagonal, delecao, insercao});

    maior.valor = 0;
    maior.salto_i = 0;
    maior.salto_j = 0;
    if (diagonal == max_score_local){
        maior.valor = diagonal;
        maior.salto_i = 1;  // salto na diagonal
        maior.salto_j = 1;
    } 

    else if (delecao == max_score_local){
        maior.valor = delecao;
        maior.salto_i = 1;  // salto para cima
        maior.salto_j = 0;
    } 

    else if (insercao == max_score_local){
        maior.valor = insercao;
        maior.salto_i = 0;
        maior.salto_j = 1; // salto da esquerda para a direita
    } 

    return maior;
}


// ----- Função que monta a matriz de alinhamento pelo método de Smith Waterman e em seguida reconstrói o alinhamento das sequências ----- //
item smith_waterman(vector<string> a, vector<string> b){
    item_sw H[(int)a.size()+1][(int)b.size()+1];
    

    // Zerando colunas especificadas
    for (int i = 0; i <= (int)a.size(); i++) {
        H[i][0].valor=0;
    }
    for (int j = 0; j <= (int)b.size(); j++) {
        H[0][j].valor = 0;
    }


    // Obtendo matriz e achando o valor máximo
    int diagonal, delecao, insercao;
    int maximo_H = 0; int max_val_i = 0; int max_val_j = 0;

    for (int i=1; i<=(int)a.size(); i++){
        for (int j=1; j<=(int)b.size(); j++){
            diagonal = H[i-1][j-1].valor + w(a[i-1],b[j-1]);
            delecao = H[i-1][j].valor - 1;
            insercao = H[i][j-1].valor - 1;

            H[i][j]=find_max(diagonal, delecao, insercao);

            if (H[i][j].valor > maximo_H) {
                    maximo_H=H[i][j].valor;
                    max_val_i=i;
                    max_val_j=j;
            }
        }
    }


// ____________________________________________ //
//  Reconstrução do alinhamento das sequências //
// __________________________________________ //

    vector<string> match_seq_a;
    vector<string> match_seq_b;

    int i=max_val_i; int j=max_val_j;


    // Reconstruindo o caminho a partir dos saltos do struct
    while ( (i>0 && j>0)  && (!(H[i][j].salto_j==0 && H[i][j].salto_i==0)) ) {
        int pos_i=i;
        int pos_j=j;
        if (H[i][j].valor == 0) break; // célula da matriz com valor zero

        if (H[pos_i][pos_j].salto_i==0 && H[pos_i][pos_j].salto_j ==1){
            match_seq_a.push_back("_");
            match_seq_b.push_back(b[j-1]);
        }
        else if (H[pos_i][pos_j].salto_i==1 && H[pos_i][pos_j].salto_j ==0){
            match_seq_a.push_back(a[i-1]);
            match_seq_b.push_back("_");
        }
        else{
            match_seq_a.push_back(a[i-1]);
            match_seq_b.push_back(b[j-1]);
        }
        
        
        i= i- H[pos_i][pos_j].salto_i;
        j=j- H[pos_i][pos_j].salto_j;

                
    }

    // Invertendo sequências
    reverse(match_seq_a.begin(),match_seq_a.end());
    reverse(match_seq_b.begin(),match_seq_b.end());


    return {maximo_H, match_seq_a, match_seq_b};
}




// Para compilar: 

// g++ -Wall -O3 -fopenmp exaustiva-paralela.cpp -o exaustiva-paralela
// ./exaustiva-paralela < dna.seq