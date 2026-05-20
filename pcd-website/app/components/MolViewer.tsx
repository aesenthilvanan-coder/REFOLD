"use client";
import { useEffect, useRef, useState } from "react";

interface MolViewerProps {
  pdbUrl: string;
  pocketResidues: string[];   // e.g. ["A40", "R41", "P42"]
  height?: number;
}

declare global {
  interface Window {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    $3Dmol: any;
  }
}

export function MolViewer({ pdbUrl, pocketResidues, height = 420 }: MolViewerProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<unknown>(null);
  const [status, setStatus] = useState<"loading" | "ready" | "error">("loading");

  useEffect(() => {
    let cancelled = false;

    async function init() {
      // Load 3Dmol.js from CDN if not already present
      if (!window.$3Dmol) {
        await new Promise<void>((resolve, reject) => {
          const s = document.createElement("script");
          s.src = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.3/3Dmol-min.js";
          s.onload = () => resolve();
          s.onerror = () => reject(new Error("3Dmol.js failed to load"));
          document.head.appendChild(s);
        });
      }

      if (cancelled || !containerRef.current) return;

      const viewer = window.$3Dmol.createViewer(containerRef.current, {
        backgroundColor: "#07090f",
        antialias: true,
      });
      viewerRef.current = viewer;

      try {
        const res = await fetch(pdbUrl);
        if (!res.ok) throw new Error(`PDB fetch failed: ${res.status}`);
        const pdbText = await res.text();
        if (cancelled) return;

        viewer.addModel(pdbText, "pdb");

        // Main protein — clean cartoon
        viewer.setStyle({}, {
          cartoon: {
            color: "#1a3a5c",
            opacity: 0.85,
            thickness: 0.3,
          }
        });

        // Pocket-lining residues — highlighted in accent cyan
        if (pocketResidues.length > 0) {
          // Parse residue numbers — use comma-separated string (most reliable format in 3Dmol.js)
          const resnums = pocketResidues
            .map(r => parseInt(r.replace(/[^0-9]/g, ""), 10))
            .filter(n => !isNaN(n) && n > 0);
          if (resnums.length > 0) {
            const resiStr = resnums.join(",");
            viewer.setStyle(
              { resi: resiStr },
              {
                cartoon: { color: "#0ea5e9", opacity: 1.0, thickness: 0.5 },
                stick: { colorscheme: "cyanCarbon", radius: 0.25, opacity: 1.0 },
              }
            );
            // Transparent surface over pocket residues
            viewer.addSurface(
              window.$3Dmol.SurfaceType.SAS,
              { opacity: 0.30, color: "#0ea5e9" },
              { resi: resiStr }
            );
          }
        }

        viewer.zoomTo();
        viewer.zoom(0.85);
        viewer.render();
        if (!cancelled) setStatus("ready");
      } catch (err) {
        console.error("MolViewer:", err);
        if (!cancelled) setStatus("error");
      }
    }

    init();
    return () => { cancelled = true; };
  }, [pdbUrl, pocketResidues]);

  return (
    <div className="relative rounded-xl overflow-hidden border border-[var(--border)] mol-viewer"
         style={{ height }}>
      <div ref={containerRef} style={{ width: "100%", height: "100%" }} />

      {/* Overlay states */}
      {status === "loading" && (
        <div className="absolute inset-0 flex items-center justify-center"
             style={{ background: "#07090f" }}>
          <div className="text-center space-y-3">
            <div className="flex justify-center gap-1">
              {[0,1,2].map(i => (
                <div key={i} className="w-1.5 h-1.5 rounded-full bg-[var(--accent)]"
                     style={{ animation: `pulse-dot 1.2s ease-in-out ${i * 0.2}s infinite` }} />
              ))}
            </div>
            <div className="section-label">Loading structure...</div>
          </div>
        </div>
      )}
      {status === "error" && (
        <div className="absolute inset-0 flex items-center justify-center"
             style={{ background: "#07090f" }}>
          <div className="text-center space-y-2">
            <div className="section-label">Structure unavailable</div>
            <p className="text-xs" style={{ color: "var(--text-muted)" }}>PDB file not found at {pdbUrl}</p>
          </div>
        </div>
      )}

      {/* Legend overlay */}
      {status === "ready" && (
        <div className="absolute bottom-3 left-3 flex items-center gap-4 px-3 py-1.5 rounded-md border border-[var(--border)]"
             style={{ background: "rgba(7,9,15,0.85)", backdropFilter: "blur(4px)" }}>
          <div className="flex items-center gap-1.5">
            <div className="w-3 h-1 rounded-sm" style={{ background: "#1a3a5c" }} />
            <span className="text-[10px]" style={{ color: "var(--text-muted)" }}>Protein</span>
          </div>
          <div className="flex items-center gap-1.5">
            <div className="w-3 h-1 rounded-sm bg-[var(--accent)]" />
            <span className="text-[10px]" style={{ color: "var(--text-muted)" }}>Pocket</span>
          </div>
        </div>
      )}
    </div>
  );
}
