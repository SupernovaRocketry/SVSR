# SR-6DOF
#######################################################################################
Contém os códigos utilizados para o simulador de vôo 6DOF da equipe Supernova Rocketry
#######################################################################################

# Códigos e Funções
O arquivo main_6DOF.m é o código principal, o qual chama todas as funções com dados de projeto e coeficientes necessários no simulador 6DOF.

O arquivo Calculos(Nakka).xls é uma planilha com todos os dados de propulsão do motor foguete classe J da equipe Supernova Rocketry.

A função dadosProjeto.m retorna todas as informações sobre o minifoguete em estruturas de dados.

A função calculaDadosAtm.m retorna os dados de densidade do ar, pressão e temperatura atmosféricas e aceleração da gravidade de acordo com a altitude.

As funções calculaAlpha.m, calculaBeta.m e calculaModulo.m retornam, respectivamente, o ângulo de ataque, o ângulo de escorregamento e o módulo de 3 componentes vetoriais.

As funções calculaCoefArrasto.m e calculaCoefSustentacao.m retornam, respectivamente, os coeficientes de arrasto e de sustentação de cada parte do minifoguete utilizados no simulador.

As funções calculaArrastoMinimo.m, calculaArrastoAoA.m e calculaCoefMomento retornam, respectivamente, o coeficiente de arrasto mínimo (desconsiderando AoA), o coeficiente de arrasto de interferência (considerando AoA) e o coeficiente de momento aerodinâmicos.

----------------------------------------------------------------------------------------------------------------------------------------
# Descrição do método de resolução

Após serem definidos todas as variáveis a serem utilizadas na função main_6DOF.m, as condições iniciais são calculadas para todas as matrizes do modelo.

Definidas as condições iniciais do sistema, é iniciado o processo de resolução do modelo determinístico do sistema dado por:

  M * dvdt + C * v + D * v + G = tau 

Para resolução da EDO que descreve o problema foi utilizado o método de Runge-Kutta de 4ª ordem.

Tendo em vista a construção do método, foi escolhida a abordagem de um loop while com restrições de suspensão do mesmo, uma vez que todos os parâmetros envolvidos são dinâmicos e, dessa forma, é possível conseguir um maior controle do laço.
