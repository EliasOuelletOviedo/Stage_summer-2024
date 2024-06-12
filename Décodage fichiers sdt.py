import io
import os
import zipfile
from typing import TYPE_CHECKING
import numpy as np
from sdtfile import SdtFile

class BlockType:
    """BLOCK_HEADER.block_type field.

    Parameters:
        value: Value of BLOCK_HEADER.block_type.

    """

    __slots__ = ('mode', 'contents', 'dtype', 'compress')

    mode: str
    """BLOCK_CREATION."""

    contents: str
    """BLOCK_CONTENT."""

    dtype: np.dtype
    """BLOCK_DTYPE."""

    compress: bool
    """Data is compressed."""

    def __init__(self, value: int, /) -> None:
        self.mode = BLOCK_CREATION[value & 0xF]
        self.contents = BLOCK_CONTENT[value & 0xF0]
        self.dtype = BLOCK_DTYPE[value & 0xF00]
        self.compress = bool(value & 0x1000)

    def __repr__(self) -> str:
        return f'<{self.__class__.__name__} {self.mode} {self.contents}>'

    def __str__(self) -> str:
        return indent(
            repr(self),
            # f'mode: {self.mode}',
            # f'contents: {self.contents}',
            f'dtype: {self.dtype}',
            f'compress: {self.compress}',
        )

def indent(*args) -> str:
    """Return joined string representations of objects with indented lines."""
    text = '\n'.join(str(arg) for arg in args)
    return '\n'.join(
        ('  ' + line if line else line) for line in text.splitlines() if line
    )[2:]


arg = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/Fichiers_sdt/4.sdt"

filename = os.fspath(arg)

fh = open(arg, "rb")

print(str(fh.read()[:1000]))

print("\n")

print(np.size(np.fromfile(fh)))

FILE_HEADER: list[tuple[str, str]] = [
    ('revision', 'i2'),
    ('info_offs', 'i4'),
    ('info_length', 'i2'),
    ('setup_offs', 'i4'),
    ('setup_length', 'u2'),
    ('data_block_offs', 'i4'),
    ('no_of_data_blocks', 'i2'),
    ('data_block_length', 'u4'),
    ('meas_desc_block_offs', 'i4'),
    ('no_of_meas_desc_blocks', 'i2'),
    ('meas_desc_block_length', 'i2'),
    ('header_valid', 'u2'),
    ('reserved1', 'u4'),
    ('reserved2', 'u2'),
    ('chksum', 'u2'),
]

SETUP_BIN_HDR: list[tuple[str, str]] = [
    ('soft_rev', 'u4'),
    ('para_length', 'u4'),
    ('reserved1', 'u4'),
    ('reserved2', 'u2'),
]

# Info collected when measurement finished
MEASURE_STOP_INFO: list[tuple[str, str]] = [
    ('status', 'u2'),
    ('flags', 'u2'),
    ('stop_time', 'f4'),
    ('cur_step', 'i4'),
    ('cur_cycle', 'i4'),
    ('cur_page', 'i4'),
    ('min_sync_rate', 'f4'),
    ('min_cfd_rate', 'f4'),
    ('min_tac_rate', 'f4'),
    ('min_adc_rate', 'f4'),
    ('max_sync_rate', 'f4'),
    ('max_cfd_rate', 'f4'),
    ('max_tac_rate', 'f4'),
    ('max_adc_rate', 'f4'),
    ('reserved1', 'i4'),
    ('reserved2', 'f4'),
]

# Info collected when FIFO measurement finished
MEASURE_FCS_INFO: list[tuple[str, str]] = [
    ('chan', 'u2'),
    ('fcs_decay_calc', 'u2'),
    ('mt_resol', 'u4'),
    ('cortime', 'f4'),
    ('calc_photons', 'u4'),
    ('fcs_points', 'i4'),
    ('end_time', 'f4'),
    ('overruns', 'u2'),
    ('fcs_type', 'u2'),
    ('cross_chan', 'u2'),
    ('mod', 'u2'),
    ('cross_mod', 'u2'),
    ('cross_mt_resol', 'u4'),
]

# Extension of MeasFCSInfo for other histograms
HIST_INFO: list[tuple[str, str]] = [
    ('fida_time', 'f4'),
    ('filda_time', 'f4'),
    ('fida_points', 'i4'),
    ('filda_points', 'i4'),
    ('mcs_time', 'f4'),
    ('mcs_points', 'i4'),
    ('cross_calc_phot', 'u4'),
    ('mcsta_points', 'u2'),
    ('mcsta_flags', 'u2'),
    ('mcsta_tpp', 'u4'),
    ('calc_markers', 'u4'),
    ('fcs_calc_phot', 'u4'),
    ('reserved3', 'u4'),
]

HIST_INFO_EXT: list[tuple[str, str]] = [
    ('first_frame_time', 'f4'),
    ('frame_time', 'f4'),
    ('line_time', 'f4'),
    ('pixel_time', 'f4'),
    ('scan_type', 'i2'),
    ('skip_2nd_line_clk', 'i2'),
    ('right_border', 'u4'),
    ('info', 'S40'),
]

MEASURE_INFO_EXT: list[tuple[str, str]] = [
    ('DCU_in_use', 'u4'),
    ('dcu_ser_no', '4S16'),
    ('axio_name', 'S32'),
    ('axio_lens_name', 'S64'),
    ('SIS_in_use', 'u4'),
    ('sis_ser_no', '4S16'),
    ('gvd_ser_no', 'S16'),
    ('gvd_zoom_factor', 'f4'),
    ('DCS_FOV_at_zoom_1', 'f4'),
    ('axio_connected', 'i2'),
    ('axio_lens_magnifier', 'f4'),
    ('axio_FOV', 'f4'),
    ('tdc_offset', '4f4'),
    ('tdc_control', 'u4'),
    ('reserve', 'S1250'),
]

# Measurement description blocks
MEASURE_INFO: list[tuple[str, str | list[tuple[str, str]]]] = [
    ('time', 'S9'),
    ('date', 'S11'),
    ('mod_ser_no', 'S16'),
    ('meas_mode', 'i2'),
    ('cfd_ll', 'f4'),
    ('cfd_lh', 'f4'),
    ('cfd_zc', 'f4'),
    ('cfd_hf', 'f4'),
    ('syn_zc', 'f4'),
    ('syn_fd', 'i2'),
    ('syn_hf', 'f4'),
    ('tac_r', 'f4'),
    ('tac_g', 'i2'),
    ('tac_of', 'f4'),
    ('tac_ll', 'f4'),
    ('tac_lh', 'f4'),
    ('adc_re', 'i2'),
    ('eal_de', 'i2'),
    ('ncx', 'i2'),
    ('ncy', 'i2'),
    ('page', 'u2'),
    ('col_t', 'f4'),
    ('rep_t', 'f4'),
    ('stopt', 'i2'),
    ('overfl', 'u1'),
    ('use_motor', 'i2'),
    ('steps', 'u2'),
    ('offset', 'f4'),
    ('dither', 'i2'),
    ('incr', 'i2'),
    ('mem_bank', 'i2'),
    ('mod', 'S16'),
    ('syn_th', 'f4'),
    ('dead_time_comp', 'i2'),
    ('polarity_l', 'i2'),
    ('polarity_f', 'i2'),
    ('polarity_p', 'i2'),
    ('linediv', 'i2'),
    ('accumulate', 'i2'),
    ('flbck_y', 'i4'),
    ('flbck_x', 'i4'),
    ('bord_u', 'i4'),
    ('bord_l', 'i4'),
    ('pix_time', 'f4'),
    ('pix_clk', 'i2'),
    ('trigger', 'i2'),
    ('scan_x', 'i4'),
    ('scan_y', 'i4'),
    ('scan_rx', 'i4'),
    ('scan_ry', 'i4'),
    ('fifo_typ', 'i2'),
    ('epx_div', 'i4'),
    ('mod_code', 'u2'),
    ('mod_fpga_ver', 'u2'),
    ('overflow_corr_factor', 'f4'),
    ('adc_zoom', 'i4'),
    ('cycles', 'i4'),
    ('StopInfo', MEASURE_STOP_INFO),
    ('FCSInfo', MEASURE_FCS_INFO),
    ('image_x', 'i4'),
    ('image_y', 'i4'),
    ('image_rx', 'i4'),
    ('image_ry', 'i4'),
    ('xy_gain', 'i2'),
    ('master_clock', 'i2'),
    ('adc_de', 'i2'),
    ('det', 'i2'),
    ('x_axis', 'i2'),
    ('MeasHISTInfo', HIST_INFO),
    ('HISTInfoExt', HIST_INFO_EXT),
    ('sync_delay', 'f4'),
    ('sdel_ser_no', 'u2'),
    ('sdel_input', 'i1'),
    ('mosaic_ctrl', 'i1'),
    ('mosaic_x', 'u1'),
    ('mosaic_y', 'u1'),
    ('frames_per_el', 'i2'),
    ('chan_per_el', 'i2'),
    ('mosaic_cycles_done', 'i4'),
    ('mla_ser_no', 'u2'),
    ('DCC_in_use', 'u1'),
    ('dcc_ser_no', 'S12'),
    ('TiSaLas_status', 'u2'),
    ('TiSaLas_wav', 'u2'),
    ('AOM_status', 'u1'),
    ('AOM_power', 'u1'),
    ('ddg_ser_no', 'S8'),
    ('prior_ser_no', 'i4'),
    ('mosaic_x_hi', 'u1'),
    ('mosaic_y_hi', 'u1'),
    ('reserve', 'S11'),
    ('extension_used', 'u1'),
    ('minfo_ext', MEASURE_INFO_EXT),
]

BLOCK_HEADER_OLD: list[tuple[str, str]] = [
    ('block_no', 'i2'),
    ('data_offs', 'i4'),
    ('next_block_offs', 'i4'),
    ('block_type', 'u2'),
    ('meas_desc_block_no', 'i2'),
    ('lblock_no', 'u4'),
    ('block_length', 'u4'),
]

BLOCK_HEADER: list[tuple[str, str]] = [
    ('data_offs_ext', 'u1'),
    ('next_block_offs_ext', 'u1'),
    ('data_offs', 'u4'),
    ('next_block_offs', 'u4'),
    ('block_type', 'u2'),
    ('meas_desc_block_no', 'i2'),
    ('lblock_no', 'u4'),
    ('block_length', 'u4'),
]

# Mode of creation
BLOCK_CREATION: dict[int, str] = {
    0: 'NOT_USED',
    1: 'MEAS_DATA',
    2: 'FLOW_DATA',
    3: 'MEAS_DATA_FROM_FILE',
    4: 'CALC_DATA',
    5: 'SIM_DATA',
    8: 'FIFO_DATA',
    9: 'FIFO_DATA_FROM_FILE',
}

BLOCK_CONTENT: dict[int, str] = {
    0x0: 'DECAY_BLOCK',
    0x10: 'PAGE_BLOCK',
    0x20: 'FCS_BLOCK',
    0x30: 'FIDA_BLOCK',
    0x40: 'FILDA_BLOCK',
    0x50: 'MCS_BLOCK',
    0x60: 'IMG_BLOCK',
    0x70: 'MCSTA_BLOCK',
    0x80: 'IMG_MCS_BLOCK',
    0x90: 'MOM_BLOCK',
    0xA0: 'IMG_INT_BLOCK',
    0xB0: 'IMG_WF_BLOCK',
    0xC0: 'IMG_LIFE_BLOCK',
}

# Data type
BLOCK_DTYPE: dict[int, np.dtype] = {
    0x000: np.dtype('<u2'),
    0x100: np.dtype('<u4'),
    0x200: np.dtype('<f8'),
}

HEADER_VALID: dict[int, bool] = {0x1111: False, 0x5555: True}

INFO_IDS: dict[str, str] = {
    'SPC Setup Script File': 'Setup script mode: setup only',
    'SPC Setup & Data File': 'Normal mode: setup + data',
    'SPC DLL Data File': 'DLL created: no setup, only data',
    'SPC Flow Data File': 'Continuous Flow mode: no setup, only data',
    'SPC FCS Data File': (
        'FIFO mode: setup, data blocks = Decay, FCS, FIDA, FILDA & MCS '
        'curves for each used routing channel'
    ),
}

with open(arg, 'rb') as fh:
    header = np.rec.fromfile(fh, dtype=FILE_HEADER, shape=1, byteorder='<')[0]
    
    fh.seek(header.data_block_offs)
    bh = np.rec.fromfile(fh, dtype=BLOCK_HEADER, shape=1, byteorder='<')[0]

    bio = io.BytesIO(fh.read(bh.next_block_offs - bh.data_offs))
    with zipfile.ZipFile(bio) as zf:
        databytes = zf.read(zf.filelist[0].filename)  # data_block
    del bio

    data = np.frombuffer(databytes, dtype="uint16", count=268435456)


print(header)
print(bh)
# print(bh.block_type)
# print(bh.block_length)
# print(len(databytes))
print(len(databytes))
# print(data[:256])

print(bh.next_block_offs - bh.data_offs)
print(bh.next_block_offs)
print(bh.data_offs)
# print(test)


# file = SdtFile(arg)

# what = []

# for i in range(100):
#     what.append(file.data[0][i][0])

# print(file.data[0][0])

# print(np.size(FILE_HEADER))