#!/bin/bash

Nodos=(1 2 4)
N=(512 1024 2048 4096)

for nodo in "${Nodos[@]}"; do
  echo "Ejecutando para ${nodo} nodos"
  for n in "${N[@]}"; do
    sbatch --tasks-per-node=8 -N "${nodo}" -o "resultados_mpi/individual-out-${n}.txt" -e "resultados_mpi/individual-error-${n}.txt" mpiBtIndividual.o "${n}";
    cat "resultados_mpi/individual-out-${n}.txt"
    echo
    # Esperando mínimo 15 segundos y pidiendo confirmación para permitir que ga$
    echo "Tiempo de gracia (15 seg)..."
    sleep 15
    echo
    read -p "Esperando input... " -n1 -s
  done;
  echo
done;
