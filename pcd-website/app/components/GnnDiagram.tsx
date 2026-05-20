"use client";

export function GnnDiagram() {
  return (
    <div className="rounded-xl overflow-hidden border border-[var(--border)]">
      {/* eslint-disable-next-line @next/next/no-img-element */}
      <img
        src="/images/arch_gnn_diffusion.png"
        alt="GNN-LM Fusion Classifier and SE(3) Diffusion Model technical architecture"
        className="w-full h-auto"
        onError={e => { e.currentTarget.style.display = "none"; }}
      />
    </div>
  );
}
