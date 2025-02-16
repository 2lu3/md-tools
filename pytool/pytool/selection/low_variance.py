from typing import Union
import MDAnalysis as mda
import warnings
import numpy as np
from loguru import logger
import japanize_matplotlib
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")


def load_position_resid(u_list: list[mda.Universe], selection: Union[str, list[str]]):
    logger.debug("Loading positions")

    if isinstance(selection, str):
        selections = [selection] * len(u_list)
    else:
        selections = selection

    positions = []
    resids = []
    for u, sel in zip(u_list, selections):
        for _ in u.trajectory:
            atoms = u.select_atoms(sel)
            positions.append(atoms.positions)
            resids.append(atoms.resids)

    positions = np.array(positions)
    return positions, resids


def kabsch(P, Q):
    """
    Kabsch アルゴリズムを用いて、点群 P を Q に最適に重ね合わせる回転 R と並進 t を計算する。
    P, Q: (N, 3) の numpy 配列
    Returns:
        R: 回転行列 (3, 3)
        t: 並進ベクトル (3,)
    """
    # 各点群の重心を計算
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)
    # 中心化
    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q
    # 共分散行列の計算
    H = np.dot(P_centered.T, Q_centered)
    # SVD 分解
    U, S, Vt = np.linalg.svd(H)
    V = Vt.T
    # 回転行列 R = V * U^T（鏡像反転を防ぐための調整）
    d = np.linalg.det(np.dot(V, U.T))
    D = np.eye(3)
    D[2, 2] = d
    R = np.dot(V, np.dot(D, U.T))
    # 並進ベクトル
    t = centroid_Q - np.dot(centroid_P, R)
    return R, t


def align_positions(positions):
    """
    全フレームの座標を、平均値に対してアライメントする。
    positions: shape (n_frames, n_residues, 3)
    Returns:
        aligned_positions: 同じ shape のアライメント後の座標
    """
    n_frames = positions.shape[0]
    aligned_positions = np.empty_like(positions)
    # 最初のフレームを参照構造とする
    ref = np.mean(positions[0], axis=0)

    for i in range(n_frames):
        P = positions[i]
        # Kabsch により P を ref に重ね合わせる変換を求める
        R, t = kabsch(P, ref)
        aligned_positions[i] = np.dot(P, R) + t
    return aligned_positions


def compute_rmsf(aligned_positions):
    """
    各残基ごとの RMSF (Root Mean Square Fluctuation) を計算する。
    aligned_positions: shape (n_frames, n_residues, 3)
    Returns:
        rmsf: (n_residues,) 各残基の RMSF 値（Å単位）
    """
    # 各残基の平均位置（フレーム平均）
    mean_positions = np.mean(aligned_positions, axis=0)  # shape: (n_residues, 3)
    # 各フレームごとの差分
    diffs = aligned_positions - mean_positions  # shape: (n_frames, n_residues, 3)
    # 各点での二乗和
    squared_diffs = np.sum(diffs**2, axis=2)  # shape: (n_frames, n_residues)
    # フレームごとの平均二乗偏差 → ルートを取る
    rmsf = np.sqrt(np.mean(squared_diffs, axis=0))  # shape: (n_residues,)
    return rmsf


def plot_delta_rmsf(num_candidates_list, avg_rmsf_list, png_path: str):
    # 平均RMSFの変化量を計算
    delta_rmsf = np.diff(avg_rmsf_list)
    # x軸は候補数の中央付近に対応させる
    x_vals = (
        np.array(num_candidates_list[:-1]) + np.array(num_candidates_list[1:])
    ) / 2

    plt.figure(figsize=(12, 8))
    plt.plot(x_vals, delta_rmsf, marker="o")
    plt.xlabel("Number of Candidate Residues")
    plt.ylabel("Delta Average RMSF (Å)")
    plt.title("RMSFの変化量")
    plt.gca().invert_xaxis()
    plt.grid(True)
    plt.savefig(png_path)
    plt.close()


def plot_rmsf(num_candidates_list, avg_rmsf_list, png_path: str):
    plt.figure(figsize=(12, 8))
    plt.plot(num_candidates_list, avg_rmsf_list, marker="o")
    plt.xlabel("Number of Candidate Residues")
    plt.ylabel("Average RMSF (Å)")
    plt.title("RMSFの候補残基数依存性")
    plt.gca().invert_xaxis()  # 候補数が減少するのでX軸を反転
    plt.grid(True)
    plt.savefig(png_path)
    plt.close()



def calc_low_variance(positions: np.ndarray) -> tuple[list[float], list[np.ndarray]]:
    candidate_mask_list: list[np.ndarray] = []
    avg_rmsf_list: list[float] = []

    # 候補残基のマスク（初期はすべて True）
    candidate_mask = np.ones(positions.shape[1], dtype=bool)
    # 候補残基が1個以上残るまで反復
    while candidate_mask.sum() > 0:
        # 現在の候補残基に対応する座標を抽出
        positions_sel = positions[
            :, candidate_mask, :
        ]  # shape: (n_frames, candidate_ount, 3)
        # アライメント：ここでは最初のフレームを参照している（必要に応じて平均構造参照への変更も可能）
        aligned_sel = align_positions(positions_sel)
        # 各候補残基ごとに RMSF を計算
        rmsf_sel = compute_rmsf(aligned_sel)  # shape: (candidate_count,)
        avg_rmsf = np.mean(rmsf_sel)

        current_candidate_count = candidate_mask.sum()

        logger.info(
            f"Candidates: {current_candidate_count}, Average RMSF: {avg_rmsf:.3f} Å"
        )

        candidate_mask_list.append(candidate_mask.copy())
        avg_rmsf_list.append(avg_rmsf)

        # 毎回、RMSF が最も高い残基を除去
        idx_to_remove = np.argmax(rmsf_sel)
        candidate_indices = np.where(candidate_mask)[0]
        remove_index = candidate_indices[idx_to_remove]
        candidate_mask[remove_index] = False

    return avg_rmsf_list, candidate_mask_list


def calc_low_variance_residues(
    u_list: list[mda.Universe], selection: Union[str, list[str]] = "protein and name CA"
):
    positions: np.ndarray
    resids: list
    positions, resids = load_position_resid(u_list, selection)

    # len(resids[i]) should be same for all i
    assert len(set(map(len, resids))) == 1, "All residues should have same length"

    avg_rmsfs, candidate_masks = calc_low_variance(positions)

    result: list[tuple[float, list[str]]] = []
    for candidate, rmsf in zip(candidate_masks, avg_rmsfs):
        residues = [resid[candidate] for resid in resids]
        result.append((rmsf, residues))

    return result


def select_low_variance_residues(
    u_list: list[mda.Universe],
    remain_res_num: int,
    selection: Union[str, list[str]] = "protein and name CA",
):
    result = calc_low_variance_residues(u_list, selection)

    residues = result[-remain_res_num][1]

    return "(resid " + " ".join(residues) + ")"
