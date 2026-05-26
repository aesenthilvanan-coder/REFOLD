import { readFile } from "fs/promises";
import path from "path";
import { PCDAtlas } from "../types";

export async function fetchAtlas(): Promise<PCDAtlas> {
  const filePath = path.join(process.cwd(), "public", "PCD_global_atlas.json");
  const raw = await readFile(filePath, "utf-8");
  return JSON.parse(raw) as PCDAtlas;
}
