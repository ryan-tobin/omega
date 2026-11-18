//! The GPU: finding one, and waiting for its answers.
//!
//! A GPU runs one small program on thousands of atoms at a time. We reach
//! it through the wgpu crate (see TECH_STACK): access to the GPU is not
//! something we write ourselves. wgpu speaks Metal on a Mac, and Vulkan or
//! DirectX elsewhere. The program the GPU runs, forces.wgsl, is ours, in
//! WGSL, the shader language of WebGPU.
//!
//! This module is built only with the `gpu` feature
//! (`cargo test --features gpu`). Without it the crate uses no third-party
//! code.
//!
//! A GPU here has no 64-bit floats: wgpu offers them on Vulkan only, not
//! on Metal. So everything it computes is in 32-bit floats and 32-bit
//! whole numbers. forces.rs says what that costs.

pub mod buffers;
pub mod forces;
pub mod sorted;

use std::future::Future;
use std::pin::pin;
use std::sync::Arc;
use std::task::{Context, Poll, Wake, Waker};
use std::thread::{self, Thread};

/// A GPU, ready for work.
pub struct Gpu {
    /// Makes buffers and programs on the GPU.
    pub device: wgpu::Device,
    /// Sends data and work to the GPU.
    pub queue: wgpu::Queue,
    /// The GPU's name, its kind, and the system library it is reached
    /// through, for the record of a benchmark.
    pub name: String,
}

impl Gpu {
    /// The GPU wgpu picks on this machine, preferring the fastest one, or
    /// None if the machine has none it can use.
    pub fn new() -> Option<Self> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
        let options = wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            ..Default::default()
        };
        let adapter = block_on(instance.request_adapter(&options)).ok()?;
        let info = adapter.get_info();
        // Ask for everything this GPU allows, such as its largest buffers.
        let descriptor = wgpu::DeviceDescriptor {
            required_limits: adapter.limits(),
            ..Default::default()
        };
        let (device, queue) = block_on(adapter.request_device(&descriptor)).ok()?;
        let name = format!("{} ({:?}, {:?})", info.name, info.device_type, info.backend);
        Some(Self {
            device,
            queue,
            name,
        })
    }
}

/// Wakes the thread that is waiting in `block_on`.
struct Unpark(Thread);

impl Wake for Unpark {
    fn wake(self: Arc<Self>) {
        self.0.unpark();
    }
}

/// Runs `future` to its end on the calling thread and returns its output.
/// wgpu hands out the GPU through futures, because in a browser the answer
/// comes later. Between polls the thread sleeps until the future wakes it.
pub fn block_on<F: Future>(future: F) -> F::Output {
    let mut future = pin!(future);
    let waker = Waker::from(Arc::new(Unpark(thread::current())));
    let mut context = Context::from_waker(&waker);
    loop {
        match future.as_mut().poll(&mut context) {
            Poll::Ready(output) => return output,
            Poll::Pending => thread::park(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::pin::Pin;
    use std::sync::atomic::{AtomicBool, Ordering};
    use std::time::Duration;

    #[test]
    fn block_on_returns_the_output_of_a_future() {
        assert_eq!(block_on(async { 6 * 7 }), 42);
    }

    /// Ready once `done` is set. The first poll starts a thread that sets
    /// it 20 ms later and then wakes the waiting thread.
    struct Later {
        done: Arc<AtomicBool>,
        started: bool,
    }

    impl Future for Later {
        type Output = u32;

        fn poll(mut self: Pin<&mut Self>, context: &mut Context<'_>) -> Poll<u32> {
            if self.done.load(Ordering::SeqCst) {
                return Poll::Ready(7);
            }
            if !self.started {
                self.started = true;
                let (done, waker) = (self.done.clone(), context.waker().clone());
                thread::spawn(move || {
                    thread::sleep(Duration::from_millis(20));
                    done.store(true, Ordering::SeqCst);
                    waker.wake();
                });
            }
            Poll::Pending
        }
    }

    #[test]
    fn block_on_sleeps_until_the_future_wakes_it() {
        let later = Later {
            done: Arc::new(AtomicBool::new(false)),
            started: false,
        };
        assert_eq!(block_on(later), 7);
    }
}
