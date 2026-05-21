"use client";
import { useEffect, useRef, useState } from "react";

interface MolViewerProps {
  pdbUrl: string;
  pocketResidues: string[];
  ligandSdfUrl?: string;
  height?: number;
}

declare global {
  interface Window {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    $3Dmol: any;
  }
}

export function MolViewer({ pdbUrl, pocketResidues, ligandSdfUrl, height = 300 }: MolViewerProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<unknown>(null);
  const [status, setStatus] = useState<"loading" | "ready" | "error">("loading");
  const [hasLigand, setHasLigand] = useState(false);

  useEffect(() => {
    let cancelled = false;

    const blockZoom = (e: WheelEvent) => e.preventDefault();
    const blockMiddle = (e: MouseEvent) => { if (e.button === 1) e.preventDefault(); };

    async function init() {
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
        backgroundColor: "#0d1424",
        antialias: true,
        disableScroll: true,
      });
      viewerRef.current = viewer;

      const el = containerRef.current;
      el.addEventListener("wheel", blockZoom, { passive: false });
      el.addEventListener("mousedown", blockMiddle);

      try {
        const res = await fetch(pdbUrl);
        if (!res.ok) throw new Error(`PDB fetch failed: ${res.status}`);
        const pdbText = await res.text();
        if (cancelled) return;

        viewer.addModel(pdbText, "pdb");

        // Protein backbone — muted navy-teal ribbon
        viewer.setStyle({}, {
          cartoon: { color: "#1e3a5c", opacity: 0.80, thickness: 0.28 }
        });

        // Pocket-lining residues — accent teal
        if (pocketResidues.length > 0) {
          const resnums: number[] = pocketResidues
            .map(r => parseInt(r.replace(/[^0-9]/g, ""), 10))
            .filter(n => !isNaN(n) && n > 0);
          if (resnums.length > 0) {
            viewer.setStyle(
              { resi: resnums },
              {
                cartoon: { color: "#0ea5e9", opacity: 1.0, thickness: 0.4 },
                stick:   { colorscheme: "cyanCarbon", radius: 0.22, opacity: 1.0 },
              }
            );
            viewer.addSurface(
              window.$3Dmol.SurfaceType.SAS,
              { opacity: 0.22, color: "#38bdf8" },
              { resi: resnums }
            );
          }
        }

        // Chaperone ligand — amber sticks
        if (ligandSdfUrl) {
          try {
            const sdfRes = await fetch(ligandSdfUrl);
            if (sdfRes.ok) {
              const sdfText = await sdfRes.text();
              if (!cancelled && sdfText.trim()) {
                viewer.addModel(sdfText, "sdf");
                const models = viewer.getModelList();
                const ligandModel = models[models.length - 1];
                ligandModel.setStyle({}, {
                  stick:  { colorscheme: "orangeCarbon", radius: 0.18, opacity: 1.0 },
                  sphere: { colorscheme: "orangeCarbon", radius: 0.26, opacity: 0.85 },
                });
                if (!cancelled) setHasLigand(true);
              }
            }
          } catch {
            // Ligand optional — protein view still valid
          }
        }

        viewer.zoomTo();
        viewer.zoom(0.82);
        viewer.render();
        if (!cancelled) setStatus("ready");
      } catch (err) {
        console.error("MolViewer:", err);
        if (!cancelled) setStatus("error");
      }
    }

    init();
    return () => {
      cancelled = true;
      const el = containerRef.current;
      if (el) {
        el.removeEventListener("wheel", blockZoom);
        el.removeEventListener("mousedown", blockMiddle);
      }
    };
  }, [pdbUrl, pocketResidues, ligandSdfUrl]);

  return (
    <div className="relative mol-viewer" style={{ height }}>
      <div ref={containerRef} style={{ width: "100%", height: "100%" }} />

      {/* Loading */}
      {status === "loading" && (
        <div className="absolute inset-0 flex items-center justify-center" style={{ background: "#0d1424" }}>
          <div className="text-center space-y-2">
            <div className="flex justify-center gap-1">
              {[0,1,2].map(i => (
                <div key={i} className="w-1 h-1 rounded-full pulse-dot"
                     style={{ background: "#0ea5e9", animationDelay: `${i * 0.2}s` }} />
              ))}
            </div>
            <div className="text-[10px] uppercase tracking-widest font-semibold" style={{ color: "#475569" }}>
              Loading structure
            </div>
          </div>
        </div>
      )}

      {/* Error */}
      {status === "error" && (
        <div className="absolute inset-0 flex items-center justify-center" style={{ background: "#0d1424" }}>
          <div className="text-center space-y-1">
            <div className="text-[10px] uppercase tracking-widest font-semibold" style={{ color: "#475569" }}>
              Structure unavailable
            </div>
          </div>
        </div>
      )}

      {/* Legend — dark semi-transparent to stay readable on dark canvas */}
      {status === "ready" && (
        <div className="absolute bottom-2 left-2 flex items-center gap-3 px-2.5 py-1 border"
             style={{
               background: "rgba(13,20,36,0.88)",
               borderColor: "rgba(255,255,255,0.08)",
               borderRadius: "2px",
             }}>
          <div className="flex items-center gap-1.5">
            <div className="w-2.5 h-1" style={{ background: "#1e3a5c", borderRadius: "1px" }} />
            <span className="text-[9px] font-semibold uppercase tracking-wider" style={{ color: "#64748b" }}>Protein</span>
          </div>
          <div className="flex items-center gap-1.5">
            <div className="w-2.5 h-1" style={{ background: "#0ea5e9", borderRadius: "1px" }} />
            <span className="text-[9px] font-semibold uppercase tracking-wider" style={{ color: "#64748b" }}>Pocket</span>
          </div>
          {hasLigand && (
            <div className="flex items-center gap-1.5">
              <div className="w-2.5 h-1" style={{ background: "#f97316", borderRadius: "1px" }} />
              <span className="text-[9px] font-semibold uppercase tracking-wider" style={{ color: "#64748b" }}>Chaperone</span>
            </div>
          )}
          <span className="text-[9px]" style={{ color: "#334155" }}>drag · rotate</span>
        </div>
      )}
    </div>
  );
}
