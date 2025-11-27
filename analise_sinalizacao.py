import csv
import sys

# Importa o algoritmo do outro arquivo
try:
    from algoritmo_dijkstra import dijkstra_sinalizacao
except ImportError:
    print("ERRO: Não foi possível importar o algoritmo do arquivo 'algoritmo_dijkstra.py'.")
    sys.exit()

def carregar_dados_string(caminho_arquivo, limite_confianca=400):
    """
    Lê o arquivo do STRING DB e monta o grafo ponderado.
    """
    grafo = {}
    print(f"Lendo base de dados: {caminho_arquivo} (com limite de confiança {limite_confianca})")
    
    try:
        with open(caminho_arquivo, 'r') as f:
            leitor = csv.reader(f, delimiter=' ')
            next(leitor) # Pula cabeçalho
            
            contagem = 0
            for linha in leitor:
                # Tratamento para linhas incompletas
                if len(linha) < 3: continue
                
                p1, p2 = linha[0], linha[1]
                try:
                    score = float(linha[2])
                except ValueError:
                    continue # Pula se o score não for número
                
                # Filtra interações fracas
                if score < limite_confianca: continue
                
                # Conversão Score -> Peso (Custo)
                # Quanto maior o score, menor o custo (caminho mais curto/fácil)
                peso = 1.0 - (score / 1000.0)
                
                if p1 not in grafo: grafo[p1] = {}
                if p2 not in grafo: grafo[p2] = {}
                
                # Armazena a aresta (grafo não direcionado)
                if p2 not in grafo[p1] or peso < grafo[p1][p2]:
                    grafo[p1][p2] = peso
                    grafo[p2][p1] = peso
                    contagem += 1
                    
        print(f"Sucesso! Grafo montado com {len(grafo)} proteínas e {contagem} conexões.")
        return grafo

    except FileNotFoundError:
        print(f"ERRO: O arquivo '{caminho_arquivo}' não foi encontrado na pasta.")
        return {}

# Lógica principal
if __name__ == "__main__":

    # Nome do arquivo baixado
    arquivo_dados = "511145.protein.links.v12.0.txt" 
    
    # Limite 700 considera apenas "Alta Confiança". 
    limite = 700 
    
    # Carrega os dados
    grafo_real = carregar_dados_string(arquivo_dados, limite_confianca=limite)
    
    if grafo_real:
        # Exibe exemplos para ajudar a escolher
        ids_disponiveis = list(grafo_real.keys())
        print("\n--- Exemplos de Proteínas disponíveis ---")
        print(f"Total: {len(ids_disponiveis)}")
        print(f"Primeiros 10 IDs: {ids_disponiveis[:10]}")
        
        # Usando os IDs dos dois primeiros da lista carregada:
        origem = ids_disponiveis[0] 
        destino = ids_disponiveis[min(10, len(ids_disponiveis)-1)] # Pega o 10º ou último
        
        print(f"\n--- Calculando Rota: {origem}  -->  {destino} ---")
        
        # Chama a função importada
        caminho, custo = dijkstra_sinalizacao(grafo_real, origem, destino)
        
        if caminho:
            print("\n" + "="*50)
            print(f"VIA DE SINALIZAÇÃO ENCONTRADA")
            print("="*50)
            print(f"Passos (Nós): {len(caminho)}")
            print(f"Custo total (Resistência): {custo:.4f}")
            print("-" * 50)
            print(" -> ".join(caminho))
            print("-" * 50)
        else:
            print("\nRESULTADO: Não existe conexão biológica conhecida entre essas proteínas")
            print("com o nível de confiança selecionado.")