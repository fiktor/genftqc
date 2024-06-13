import ctypes
import os
import sys
import functools
import numpy as np

@functools.cache
def _get_lib():
    extension = {
        "win32": ".dll",
        "darwin": ".dylib"
    }.get(sys.platform, ".so")
    lib_path = os.path.join(
        os.path.dirname(__file__),
        'lib',
        'libmultiply_gates' + extension)
    lib = ctypes.CDLL(lib_path)
    return lib

@functools.cache
def _get_cayley_lib():
    lib = _get_lib()

    lib.create_cayley_graph.argtypes = []
    lib.create_cayley_graph.restype = ctypes.c_void_p

    lib.destroy_cayley_graph.argtypes = [ctypes.c_void_p]
    lib.destroy_cayley_graph.restype = None

    lib.cayley_graph_read.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
    lib.cayley_graph_read.restype = ctypes.c_bool

    lib.cayley_graph_fill_downstream_ids.argtypes = [ctypes.c_void_p]
    lib.cayley_graph_fill_downstream_ids.restype = None

    lib.create_cycle_tops.argtypes = [ctypes.c_void_p]
    lib.create_cycle_tops.restype = ctypes.c_void_p

    lib.destroy_cycle_tops.argtypes = [ctypes.c_void_p]
    lib.destroy_cycle_tops.restype = None

    lib.cycle_tops_get_num_vertices.argtypes = [ctypes.c_void_p]
    lib.cycle_tops_get_num_vertices.restype = ctypes.c_int32

    lib.cycle_tops_get_vertices.argtypes = [
            ctypes.c_void_p, ctypes.POINTER(ctypes.c_uint32)]
    lib.cycle_tops_get_vertices.restype = None

    lib.cycle_tops_get_num_edges.argtypes = [ctypes.c_void_p]
    lib.cycle_tops_get_num_edges.restype = ctypes.c_int32

    lib.cycle_tops_get_edges.argtypes = [
            ctypes.c_void_p, ctypes.POINTER(ctypes.c_uint32),
            ctypes.POINTER(ctypes.c_uint32), ctypes.POINTER(ctypes.c_uint8)]
    lib.cycle_tops_get_edges.restype = None

    lib.create_cycle_builder.argtypes = [ctypes.c_void_p, ctypes.c_uint64]
    lib.create_cycle_builder.restype = ctypes.c_void_p

    lib.destroy_cycle_builder.argtypes = [ctypes.c_void_p]
    lib.destroy_cycle_builder.restype = None

    lib.cycle_builder_get_max_cycle_size.argtypes = [ctypes.c_void_p]
    lib.cycle_builder_get_max_cycle_size.restype = ctypes.c_uint32

    lib.cycle_builder_get_cycle_gates.argtypes = [
            ctypes.c_void_p, ctypes.c_bool,
            ctypes.POINTER(ctypes.c_uint8), ctypes.c_int32]
    lib.cycle_builder_get_cycle_gates.restype = ctypes.c_int32

    lib.cycle_builder_get_cycle_vertices.argtypes = [
            ctypes.c_void_p, ctypes.POINTER(ctypes.c_uint32), ctypes.c_int32]
    lib.cycle_builder_get_cycle_vertices.restype = ctypes.c_int32

    lib.cycle_builder_from_vertex.argtypes = [ctypes.c_void_p, ctypes.c_uint32]
    lib.cycle_builder_from_vertex.restype = None

    lib.cycle_builder_from_edge.argtypes = [
            ctypes.c_void_p, ctypes.c_uint32, ctypes.c_uint32, ctypes.c_uint8]
    lib.cycle_builder_from_edge.restype = None

    lib.cycle_builder_from_random_cycle_top.argtypes = [
            ctypes.c_void_p, ctypes.c_void_p]
    lib.cycle_builder_from_random_cycle_top.restype = None

    lib.create_openqasm_generator.argtypes = [ctypes.c_void_p]
    lib.create_openqasm_generator.restype = ctypes.c_void_p

    lib.destroy_openqasm_generator.argtypes = [ctypes.c_void_p]
    lib.destroy_openqasm_generator.restype = None

    lib.cycle_to_openqasm.argtypes = [
        ctypes.c_void_p, ctypes.POINTER(ctypes.c_uint8),
        ctypes.c_int32, ctypes.c_char_p, ctypes.c_int32]
    lib.cycle_to_openqasm.restype = ctypes.c_int32

    lib.cycle_to_tokens.argtypes = [
            ctypes.c_void_p, ctypes.POINTER(ctypes.c_uint8), ctypes.c_int32,
            ctypes.POINTER(ctypes.c_uint32), ctypes.c_int32]
    lib.cycle_to_tokens.restype = ctypes.c_int32

    lib.detokenize.argtypes = [
        ctypes.c_void_p, ctypes.POINTER(ctypes.c_uint32), ctypes.c_int32,
        ctypes.POINTER(ctypes.c_char), ctypes.c_int32]
    lib.detokenize.restype = ctypes.c_int32

    lib.gen_multiply_gates_tasks.argtypes = [
        # OpenQasmGenerator *, CycleBuilderWrapper *, CycleTops *:
        ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
        # split, num_tasks, max_size:
        ctypes.c_int32, ctypes.c_int32, ctypes.c_int32,
        # part_lengths, task_tokens:
        ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_uint32)]
    # Returns 0 on success:
    lib.gen_multiply_gates_tasks.restype = ctypes.c_int32

    return lib

lib = _get_cayley_lib()

# Must be consistent with SplitType enum in cayley_graph.h:
SPLIT_ALL = 0
SPLIT_TRAIN = 1
SPLIT_TEST = 2

class CayleyGraph:
    def __init__(self):
        self.obj = lib.create_cayley_graph()
        self._gen = None
        if not self.obj:
            raise RuntimeError("Failed to create CayleyGraph object")

    def __del__(self):
        if self.obj:
            lib.destroy_cayley_graph(self.obj)
        if self._gen:
            lib.destroy_openqasm_generator(self._gen)

    def read(self, filename):
        result = lib.cayley_graph_read(self.obj, filename.encode('utf-8'))
        return bool(result)

    def fill_downstream_ids(self):
        lib.cayley_graph_fill_downstream_ids(self.obj)

    def _get_gen(self):
        """
        Get the OpenQasmGenerator object.

        Note: self._gen should not be used directly outside __init__, __del__,
        _get_gen because it is not guaranteed to be initialized.
        """
        if not self._gen:
            self._gen = lib.create_openqasm_generator(self.obj)
            if not self._gen:
                raise RuntimeError("Failed to create OpenQasmGenerator object")
        return self._gen

    def to_openqasm(self, gates, out_size=0):
        gen = self._get_gen()
        if gates.dtype != np.uint8:
            raise ValueError("gates must be of dtype np.uint8")

        gates_ptr = gates.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8))
        if out_size <= 0:
            out_size = 18 * (len(gates) + 2)
        out_buffer = ctypes.create_string_buffer(out_size)
        result_size = lib.cycle_to_openqasm(
                gen, gates_ptr, len(gates), out_buffer, out_size)
        if result_size == -1:
            # TODO: return truncated output and a flag.
            raise RuntimeError("Output buffer size was too small")
        return out_buffer.value.decode('ascii')

    def to_tokens(self, gates, out_size=0):
        gen = self._get_gen()
        if gates.dtype != np.uint8:
            raise ValueError("gates must be of dtype np.uint8")

        gates_ptr = gates.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8))
        if out_size <= 0:
            # A 3-qubit gate has 4 tokens: token name and 3 qubit indices.
            out_size = 4 * (len(gates) + 1)
        out_buffer = np.empty(out_size, dtype=np.uint32)
        result_size = lib.cycle_to_tokens(
                gen, gates_ptr, len(gates),
                out_buffer.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                out_size)
        if result_size == -1:
            # TODO: return truncated output and a flag.
            raise RuntimeError("Output buffer size was too small")
        return out_buffer[:result_size]

    def detokenize(self, tokens, out_size=0):
        gen = self._get_gen()
        if tokens.dtype != np.uint32:
            raise ValueError("tokens must be of dtype np.uint32")

        tokens_ptr = tokens.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
        if out_size <= 0:
            out_size = 18 * (len(tokens) + 1)
        out_buffer = ctypes.create_string_buffer(out_size)
        result_size = lib.detokenize(
                gen, tokens_ptr, len(tokens), out_buffer, out_size)
        if result_size == -1:
            # TODO: return truncated output and a flag.
            raise RuntimeError("Output buffer size was too small")
        return out_buffer.value.decode('ascii')

    def gen_multiply_gates_tasks(
            self, builder, cycle_tops,
            split=SPLIT_TRAIN, num_tasks=1024, max_size=43):
        """
        Generate multiply gates tasks.

        The tasks are represented as 2 arrays:
        - part_lengths: (num_tasks, 2) array of dtype np.int32.
            part_lengths[j, 0] is the length of the "question" part of task j.
            part_lengths[j, 1] is the length of the "answer" part of task j.
        - task_tokens: (num_tasks, max_size) array of dtype np.uint32.
        """
        assert num_tasks > 0
        part_lengths = np.empty(num_tasks * 2, dtype=np.int32)
        task_tokens = np.empty(num_tasks * max_size, dtype=np.uint32)
        gen = self._get_gen()
        # gen_multiply_gates_tasks fills the unused space of task_tokens
        # with 0s (= kTokenIdEos).
        ret = lib.gen_multiply_gates_tasks(
            gen, builder.obj, cycle_tops.obj, split, num_tasks, max_size,
            part_lengths.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
            task_tokens.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)))
        assert ret == 0, f"Error generating tasks: {ret}"
        return (
            part_lengths.reshape((num_tasks, 2)),
            task_tokens.reshape((num_tasks, max_size)))

class CycleTops:
    def __init__(self, cayley_graph):
        self.obj = lib.create_cycle_tops(cayley_graph.obj)
        if not self.obj:
            raise RuntimeError("Failed to create CycleTops object")

    def __del__(self):
        if self.obj:
            lib.destroy_cycle_tops(self.obj)

    def get_vertices(self):
        num_vertices = lib.cycle_tops_get_num_vertices(self.obj)
        buffer = np.empty(num_vertices, dtype=np.uint32)
        lib.cycle_tops_get_vertices(
            self.obj, buffer.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)))
        return buffer

    def get_edges(self):
        num_edges = lib.cycle_tops_get_num_edges(self.obj)
        buffer_u = np.empty(num_edges, dtype=np.uint32)
        buffer_v = np.empty(num_edges, dtype=np.uint32)
        buffer_g = np.empty(num_edges, dtype=np.uint8)
        lib.cycle_tops_get_edges(
            self.obj,
            buffer_u.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
            buffer_v.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
            buffer_g.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)))
        return buffer_u, buffer_v, buffer_g

class CycleBuilder:
    def __init__(self, cayley_graph, seed):
        self.obj = lib.create_cycle_builder(cayley_graph.obj, seed)
        if not self.obj:
            raise RuntimeError("Failed to create CycleBuilder object")

    def __del__(self):
        if self.obj:
            lib.destroy_cycle_builder(self.obj)

    def get_max_cycle_size(self):
        return lib.cycle_builder_get_max_cycle_size(self.obj)

    def get_cycle_gates(self, invert_last=False, max_size=None):
        if max_size is None:
            max_size = self.get_max_cycle_size()
        buffer = np.empty(max_size, dtype=np.uint8)
        actual_size = lib.cycle_builder_get_cycle_gates(
            self.obj,
            ctypes.c_bool(invert_last),
            buffer.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)),
            max_size)
        if actual_size < 0:
            raise RuntimeError("Error getting cycle gates")
        return buffer[:actual_size]

    def get_cycle_vertices(self, max_size):
        buffer = (ctypes.c_uint32 * max_size)()
        actual_size = lib.cycle_builder_get_cycle_vertices(self.obj, buffer, max_size)
        if actual_size < 0:
            raise RuntimeError("Error getting cycle vertices")
        return buffer[:actual_size]

    def from_vertex(self, vertex_id):
        lib.cycle_builder_from_vertex(self.obj, vertex_id)

    def from_edge(self, u, v, g):
        lib.cycle_builder_from_edge(self.obj, u, v, g)

    def from_random_cycle_top(self, cycle_tops):
        lib.cycle_builder_from_random_cycle_top(self.obj, cycle_tops.obj)
