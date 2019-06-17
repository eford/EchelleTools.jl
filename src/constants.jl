# Physical constants
const speed_of_light_mps = 299792458.0 # m/s
export speed_of_light_mps

# Limit data size/memory usage
const max_num_pixels_per_order = 8096
const max_num_orders = 128
const max_template_length = max_num_orders * max_num_pixels_per_order
const max_num_obs = 1024
const max_num_basis_vectors = 16
const min_num_star_basis_vectors = 1
const max_num_star_basis_vectors = max_num_basis_vectors
const min_num_telluric_basis_vectors = 0
const max_num_telluric_basis_vectors = max_num_basis_vectors
export max_num_pixels_per_order, max_num_orders, max_template_length, max_num_obs
export max_num_basis_vectors, min_num_star_basis_vectors, max_num_star_basis_vectors, min_num_telluric_basis_vectors, max_num_telluric_basis_vectors


# TODO: move to each instrument?
const default_blaze_degree = 4
const default_pixel_buffer_lo = 0
const default_pixel_buffer_hi = 0
export default_blaze_degree, default_pixel_buffer_lo, default_pixel_buffer_hi
