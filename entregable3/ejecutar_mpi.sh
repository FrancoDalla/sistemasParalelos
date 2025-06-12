#!/bin/bash

Nodos=(2 4)
N=(512 1024 2048 4096)

for nodo in "${Nodos[@]}"; do
  echo "Ejecutando para ${nodo} nodos"
  for n in "${N[@]}"; do
    sbatch --tasks-per-node=1 -N "${nodo}" --wait -o "hibrido/out-2-${n}.txt" -$
    cat "hibrido/out-2-${n}.txt"
    echo
    # Esperando mínimo 15 segundos y pidiendo confirmación para permitir que ga$
    echo "Tiempo de gracia (15 seg)..."
    sleep 15
    echo
    read -p "Esperando input... " -n1 -s
  done;
  echo
done;

