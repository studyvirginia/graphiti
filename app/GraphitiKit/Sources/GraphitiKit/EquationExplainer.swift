import Foundation
import Security

public enum EquationExplainer {
    static let model = "claude-haiku-5-5"

    public enum Failure: Error, LocalizedError {
        case noKey, http(Int, String), empty
        public var errorDescription: String? {
            switch self {
            case .noKey: "Add your Anthropic API key in You › Settings to get an explanation."
            case .http(let code, let msg): "Claude API error \(code): \(msg)"
            case .empty: "Claude returned no text."
            }
        }
    }

    public static func explain(latex: String, stats: SolveSession.Stats) async throws -> String {
        guard let key = APIKeyStore.key, !key.isEmpty else { throw Failure.noKey }
        let prompt = """
        A person drew a curve with their finger and a symbolic-regression engine found this equation for it \
        (LaTeX, constants rounded to one decimal): \(latex)
        Search stats: \(String(format: "%.1f", stats.seconds)) s, \(stats.generations) shown generations, \
        \(stats.equationsTried) candidate equations tried, \(stats.terms) additive terms, shape fit \(String(format: "%.1f%%", stats.fit * 100)).
        In under 120 words of plain text (no Markdown, no LaTeX), say what each part of the equation does to the shape \
        (if it is a parametric pair x(t), y(t), describe the path it traces), then one sentence putting the stats in \
        context. Speak to someone who is not a mathematician.
        """
        var req = URLRequest(url: URL(string: "https://api.anthropic.com/v1/messages")!)
        req.httpMethod = "POST"
        req.timeoutInterval = 30
        req.setValue(key, forHTTPHeaderField: "x-api-key")
        req.setValue("2023-06-01", forHTTPHeaderField: "anthropic-version")
        req.setValue("application/json", forHTTPHeaderField: "content-type")
        req.httpBody = try JSONSerialization.data(withJSONObject: [
            "model": model, "max_tokens": 600,
            "messages": [["role": "user", "content": prompt]],
        ])
        let (data, resp) = try await URLSession.shared.data(for: req)
        let json = (try? JSONSerialization.jsonObject(with: data)) as? [String: Any]
        let code = (resp as? HTTPURLResponse)?.statusCode ?? 0
        guard code == 200 else {
            let msg = ((json?["error"] as? [String: Any])?["message"] as? String) ?? HTTPURLResponse.localizedString(forStatusCode: code)
            Diag.log("explain", "HTTP \(code): \(msg)", level: .error)
            throw Failure.http(code, msg)
        }
        let blocks = json?["content"] as? [[String: Any]] ?? []
        let text = blocks.filter { $0["type"] as? String == "text" }.compactMap { $0["text"] as? String }.joined()
        guard !text.isEmpty else { throw Failure.empty }
        return text.trimmingCharacters(in: .whitespacesAndNewlines)
    }
}

public enum APIKeyStore {
    private static let base: [String: Any] = [kSecClass as String: kSecClassGenericPassword,
                                              kSecAttrService as String: "com.graphart.anthropic",
                                              kSecAttrAccount as String: "api-key"]
    public static var key: String? {
        var q = base; q[kSecReturnData as String] = true; q[kSecMatchLimit as String] = kSecMatchLimitOne
        var out: CFTypeRef?
        guard SecItemCopyMatching(q as CFDictionary, &out) == errSecSuccess, let d = out as? Data else { return nil }
        return String(data: d, encoding: .utf8)
    }
    public static func set(_ value: String?) {
        SecItemDelete(base as CFDictionary)
        guard let v = value?.trimmingCharacters(in: .whitespacesAndNewlines), !v.isEmpty else { return }
        var q = base; q[kSecValueData as String] = Data(v.utf8)
        q[kSecAttrAccessible as String] = kSecAttrAccessibleWhenUnlockedThisDeviceOnly
        SecItemAdd(q as CFDictionary, nil)
    }
}
