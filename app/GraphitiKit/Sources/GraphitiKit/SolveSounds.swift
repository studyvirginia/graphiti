import AVFoundation

@MainActor
public final class SolveSounds {
    public static let shared = SolveSounds()

    public enum Cue: Hashable { case term(Int), tune(Int), prune, lock }

    private let engine = AVAudioEngine()
    private let player = AVAudioPlayerNode()
    private let format = AVAudioFormat(standardFormatWithSampleRate: 44_100, channels: 1)!
    private var started = false
    private var cache: [Cue: AVAudioPCMBuffer] = [:]

    private let scale: [Double] = [523.25, 587.33, 659.25, 783.99, 880.00, 1046.50, 1174.66, 1318.51, 1567.98]

    public var enabled: Bool { UserDefaults.standard.object(forKey: "soundsOn") as? Bool ?? true }

    public func play(_ cue: Cue) {
        guard enabled else { return }
        if !started { start() }
        guard started, let buffer = buffer(for: cue) else { return }
        player.scheduleBuffer(buffer, at: nil, options: [], completionHandler: nil)
    }

    private func start() {
        #if os(iOS)
        try? AVAudioSession.sharedInstance().setCategory(.ambient, options: [.mixWithOthers])
        try? AVAudioSession.sharedInstance().setActive(true)
        #endif
        engine.attach(player)
        engine.connect(player, to: engine.mainMixerNode, format: format)
        engine.mainMixerNode.outputVolume = 0.6
        do { try engine.start(); player.play(); started = true } catch { started = false }
    }

    private func buffer(for cue: Cue) -> AVAudioPCMBuffer? {
        if let b = cache[cue] { return b }
        let b: AVAudioPCMBuffer?
        switch cue {
        case .term(let i):
            b = tone([(scale[min(i, scale.count - 1)], 1), (scale[min(i, scale.count - 1)] * 2, 0.25)],
                     duration: 0.18, decay: 22, gain: 0.16)
        case .tune(let i):
            b = tone([(scale[min(i, scale.count - 1)] * 2, 1)], duration: 0.06, decay: 60, gain: 0.07)
        case .prune:
            b = tone([(196, 1), (147, 0.6)], duration: 0.16, decay: 26, gain: 0.14)
        case .lock:
            b = tone([(523.25, 1), (659.25, 0.8), (783.99, 0.7), (1046.5, 0.5), (2093, 0.12)],
                     duration: 1.1, decay: 4.2, gain: 0.11)
        }
        cache[cue] = b
        return b
    }

    private func tone(_ partials: [(Double, Double)], duration: Double, decay: Double, gain: Double) -> AVAudioPCMBuffer? {
        let frames = AVAudioFrameCount(duration * format.sampleRate)
        guard let buffer = AVAudioPCMBuffer(pcmFormat: format, frameCapacity: frames),
              let data = buffer.floatChannelData?[0] else { return nil }
        buffer.frameLength = frames
        let total = partials.reduce(0) { $0 + $1.1 }
        for n in 0..<Int(frames) {
            let t = Double(n) / format.sampleRate
            let envelope = min(1, t / 0.004) * exp(-t * decay)
            var s = 0.0
            for (f, a) in partials { s += sin(2 * .pi * f * t) * a }
            data[n] = Float(s / total * envelope * gain)
        }
        return buffer
    }
}
