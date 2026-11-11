"""Resource guard only; no mathematical/simulation code is shared by engines."""
import ctypes
import atexit,json,os,sys,time
KEEP=[]
def guard():
    if os.name!='nt':
        return {'pid':os.getpid(),'script':os.path.basename(sys.argv[0]),'memory_guard':'not applied (not Windows)'}
    from ctypes import wintypes
    class MEM(ctypes.Structure):
        _fields_=[('length',wintypes.DWORD),('load',wintypes.DWORD)]+[(x,ctypes.c_ulonglong) for x in ('total','avail','page','availpage','virt','availvirt','extended')]
    m=MEM(); m.length=ctypes.sizeof(m)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m)): raise OSError('memory query')
    if m.avail < 9*1024**3: raise RuntimeError('less than 9 GiB free')
    class BASIC(ctypes.Structure):
        _fields_=[('process_time',ctypes.c_longlong),('job_time',ctypes.c_longlong),('flags',wintypes.DWORD),('minws',ctypes.c_size_t),('maxws',ctypes.c_size_t),('active',wintypes.DWORD),('affinity',ctypes.c_size_t),('priority',wintypes.DWORD),('scheduling',wintypes.DWORD)]
    class IO(ctypes.Structure):
        _fields_=[(x,ctypes.c_ulonglong) for x in ('read_ops','write_ops','other_ops','read_bytes','write_bytes','other_bytes')]
    class EXT(ctypes.Structure):
        _fields_=[('basic',BASIC),('io',IO),('process_memory',ctypes.c_size_t),('job_memory',ctypes.c_size_t),('peak_process',ctypes.c_size_t),('peak_job',ctypes.c_size_t)]
    kernel=ctypes.windll.kernel32
    kernel.CreateJobObjectW.restype=wintypes.HANDLE
    kernel.GetCurrentProcess.restype=wintypes.HANDLE
    kernel.SetInformationJobObject.argtypes=[wintypes.HANDLE,ctypes.c_int,ctypes.c_void_p,wintypes.DWORD]
    kernel.AssignProcessToJobObject.argtypes=[wintypes.HANDLE,wintypes.HANDLE]
    job=kernel.CreateJobObjectW(None,None)
    ext=EXT(); ext.basic.flags=0x2300; ext.process_memory=1000000000;ext.job_memory=1000000000
    if not kernel.SetInformationJobObject(job,9,ctypes.byref(ext),ctypes.sizeof(ext)): raise ctypes.WinError()
    if not kernel.AssignProcessToJobObject(job,kernel.GetCurrentProcess()): raise ctypes.WinError()
    KEEP.append(job)
    record={'pid':os.getpid(),'script':os.path.basename(sys.argv[0]),'free_bytes_before':m.avail,'process_cap_bytes':1000000000,'job_cap_bytes':1000000000}
    start=time.monotonic()
    def finish():
        current=EXT()
        kernel.QueryInformationJobObject.argtypes=[wintypes.HANDLE,ctypes.c_int,ctypes.c_void_p,wintypes.DWORD,ctypes.c_void_p]
        if not kernel.QueryInformationJobObject(job,9,ctypes.byref(current),ctypes.sizeof(current),None):raise ctypes.WinError()
        record.update(peak_process_bytes=current.peak_process,peak_job_bytes=current.peak_job,wall_seconds=round(time.monotonic()-start,6))
        path=os.path.join(os.path.dirname(os.path.abspath(__file__)),os.path.splitext(os.path.basename(sys.argv[0]))[0]+('_'+'_'.join(os.path.splitext(os.path.basename(arg))[0] for arg in sys.argv[1:]) if len(sys.argv)>1 else '')+'.run.json')
        with open(path,'w') as f:json.dump(record,f,indent=2)
    atexit.register(finish)
    return record




